from collections import defaultdict
from fractions import Fraction
from functools import reduce
from itertools import product
from operator import mul

from constants import *
from exceptions import InvalidState, NonMendelianPattern
from punnet_square import punnet_square


def affected_genotype(genotype, mode):
    """Determine if the genotype will have an affected phenotype

    Args:
        genotype: Eg. 'Aa' / 'Xy'
        mode: mode of inheritance
        gender: 'MALE' or 'FEMALE'

    Returns: bool
    """
    if mode == AUTOSOMAL_DOMINANT:
        if 'A' in genotype:
            return True
    elif mode == AUTOSOMAL_RECESSIVE:
        if 'aa' in genotype:
            return True
    elif mode == X_LINKED_DOMINANT:
        if 'X' in genotype:
            return True
    elif mode == X_LINKED_RECESSIVE and 'y' in genotype.lower():
        if 'x' in genotype:
            return True
    elif mode == X_LINKED_RECESSIVE and 'y' not in genotype.lower():
        if 'xx' in genotype:
            return True

    elif mode == Y_LINKED and 'y' in genotype.lower():
        return True
    return False


def phenotype_given_genotype(genotype, mode):
    """Return the possible phenotypes for a given genotype and mode of inheritance
    """
    if mode == AUTOSOMAL_DOMINANT:
        if 'A' in genotype:
            return 'FM'
        else:
            return 'fm'
    elif mode == AUTOSOMAL_RECESSIVE:
        if 'aa' in genotype:
            return 'FM'
        else:
            return 'fm'
    elif mode == X_LINKED_DOMINANT:
        if 'y' in genotype.lower():
            if 'X' in genotype:
                return 'M'
            else:
                return 'm'
        else:
            if 'X' in genotype:
                return 'F'
            else:
                return 'f'
    elif mode == X_LINKED_RECESSIVE:
        if 'y' in genotype.lower():
            if 'x' in genotype:
                return 'M'
            else:
                return 'm'
        else:
            if 'xx' in genotype:
                return 'F'
            else:
                return 'f'

    elif mode == Y_LINKED:
        if 'y' in genotype.lower():
            return 'M'
        else:
            return 'f'


def genotype_possibilities(mode, gender, affected=None, observation=None):
    """Return valid genotypes for a given observation and mode

    This is equivalent to a null prior.

    Args:
        mode: The mode of inheritance
        gender: Optionally supply the gender
        affected: True the genotypes will be affected for the given mode

    Returns: list of valid genotypes
    """

    if mode == Y_LINKED and gender == FEMALE and affected == True:
        raise InvalidState('Can not have an affected female if condition is Y_LINKED')

    if observation is not None:
        affected = observation.affected

    if mode in AUTOSOMAL_MODES:
        chromosomes = AUTOSOMAL_CHROMOSOMES
    elif gender:
        chromosomes = {MALE: MALE_SEX_LINKED_CHROMOSOMES,
                       FEMALE: FEMALE_SEX_LINKED_CHROMOSOMES}[gender]
    else:
        chromosomes = MALE_SEX_LINKED_CHROMOSOMES.union(FEMALE_SEX_LINKED_CHROMOSOMES)

    if gender is not None and affected is not None:
        chromosomes = {genotype for genotype in chromosomes if
                       affected_genotype(genotype, mode) == affected}

    prior = Fraction(1, len(chromosomes))

    return {chromosome: prior for chromosome in chromosomes}


def normalize_probabilities(probabilities):
    """Normalize a set of probabilities.

    This is equivalent to applying bayes rule to all probabilities
    with a null prior.

    Args:
        probabilities: dict of probabilities to normalize

    Returns: dict of probabilities
    """
    prior = Fraction(1, len(probabilities))
    normalize = sum([i * prior for i in probabilities.values()])
    return {key: (value * prior) / normalize for key, value in probabilities.items()
            if value != 0}


def constrain_probabilities(mode, observation, genotype_probabilities):
    """Given a dict of genotype probabilities remove all genotypes that
    don't match the observation given a mode of inheritance.

    Args:
        mode: The mode of inheritance
        observation: Observation to match
        genotype_probabilities: dict of genotype probabilities

    Returns: dict of probabilities
    """
    constrained = {genotype: probability for genotype, probability in genotype_probabilities.items()
                   if observation in phenotype_given_genotype(genotype, mode)}
    if constrained:
        constrained = normalize_probabilities(constrained)
    return constrained


def punnet_occurrences(mode, observation, genotype_probabilities):
    """Calculate the times the observation would have been made given a dict of
    genotype probabilities.

    genotype_probabilities will generally be the result of a punnet square.

    Args:
        mode: The mode of inheritance
        observation: The observation to match
        genotype_probabilities: A dict of genotype probabilities

    Returns: Fraction indicating probability
    """
    return sum(probability for genotype, probability in genotype_probabilities.items()
               if observation in phenotype_given_genotype(genotype, mode))


def _next_parent_probabilities(mode, child, genotype_probabilities):
    """Create the next mother and father genotypes for the parent probabilities
    function"""
    mother, father = child.mother_father
    constrained_probabilities = constrain_probabilities(mode, child, genotype_probabilities)
    mother_genotypes = father_genotypes = None
    if mother is child:
        mother_genotypes = constrained_probabilities
        father_genotypes = genotype_possibilities(mode, MALE, observation=father)
    if father is child:
        father_genotypes = constrained_probabilities
        mother_genotypes = genotype_possibilities(mode, FEMALE, observation=mother)
    return mother_genotypes, father_genotypes


def parent_likelihood_n_children(mode, children, mother_genotypes, father_genotypes):
    """Return a probability that reflects how well the mother and father genotypes
    explain the observed children.

    mother and father genotypes are dicts of genotype probabilities. This is because
    we can iterate through each individual parent genotype at the top level, but any
    children will have a range of probabilities that needs to be allowed for with the
    recursive call to self

    Args:
        mode: The mode of inheritance
        observation: The observation to start from
        mother_genotypes: The probabilities of mother genotypes
        father_genotypes: the probabilities of father genotypes

    """
    genotypes_probabilities = punnet_square(mother_genotypes, father_genotypes)
    result = [1]
    for child in children:
        result.append(punnet_occurrences(mode, child, genotypes_probabilities))
        result.append(
            parent_likelihood_n_children(mode, child.children,
                                         *_next_parent_probabilities(mode, child, genotypes_probabilities))
        )
    return reduce(mul, result)


def parent_probabilities_n_children(mode, children, mother_genotypes, father_genotypes):
    """Calculate the probabilities for each possible parent combination from the supplied
    mother and father genotype probabilities

    Args
        mode: The mode of inheritance
        children: An iterable of children observations
        mother_genotype: dict of mother genotype probabilities
        father_genotype: dict of father genotype probabilities

    Returns: dict of parent combination probabilities keys in the form of (mother, father) genotype
    """
    return normalize_probabilities({
        (m_genotype, f_genotype):
            parent_likelihood_n_children(mode, children, m_genotype, f_genotype) * (f_prob * m_prob)
        for (m_genotype, m_prob), (f_genotype, f_prob) in
        product(mother_genotypes.items(), father_genotypes.items())
    })


def parent_to_children_probabilities(parent_probabilities):
    """Convert the probabilities of parent combinations to there resulting
    children probabilities.

    Args:
        parent_probabilities: dict of probabilities with keys in the form of (mother, father) genotype

    Returns: dict of genotype probabilities
    """
    child_probabilities = defaultdict(lambda: 0)
    for (f_genotype, m_genotype), parent_probability in parent_probabilities.items():
        for child, child_probability in punnet_square(f_genotype, m_genotype).items():
            child_probabilities[child] += child_probability * parent_probability

    return normalize_probabilities(child_probabilities)


def split_parent_probabilities(probabilities):
    mother = defaultdict(lambda: 0)
    father = defaultdict(lambda: 0)
    for (m_genotype, f_genotype), probability in probabilities.items():
        mother[m_genotype] += probability
        father[f_genotype] += probability
    return dict(mother), dict(father)


def back_propagate(mode, observation, mother_genotypes, father_genotypes):
    """Back propagate the information of the downstream observations to the current observation. Eg.

    What is the current observation genotype probabilities given how well the observations parent
    probabilities match the parents children?

    Args:
        mode: Mode of inheritance
        observation: The observation in question
        mother_genotypes: A dict of genotype probabilities for the observations mother
        father_genotypes: A dict of genotype probabilities for the observations father

    Returns: A dict of genotype probabilities for the observation
    """
    a = parent_probabilities_n_children(mode, observation.siblings, mother_genotypes, father_genotypes)
    b = parent_to_children_probabilities(a)
    return constrain_probabilities(mode, observation, b)


def parent_probabilities(mode, mother, father):
    """Calculate the best fit genotype probabilities for the supplied mother and father
    observation

    Args:
        mode: The mode of inheritance
        mother: The observed mother
        father: The observed father

    Returns: tuple of (mother_genotypes, father_genotypes)
    """
    mm = mf = fm = ff = None

    if mother:
        mm = mother.mother
        mf = mother.father
    if father:
        fm = father.mother
        ff = father.father

    mother_genotypes, father_genotypes = genotype_possibilities(mode, FEMALE, observation=mother), \
                                         genotype_possibilities(mode, MALE, observation=father)

    # This is how we account for missing parents.
    if mother is None and father is None:
        return mother_genotypes, father_genotypes
    elif mother is None:
        return mother_genotypes, back_propagate(mode, father, *parent_probabilities(mode, fm, ff))
    elif father is None:
        return back_propagate(mode, mother, *parent_probabilities(mode, mm, mf)), father_genotypes
    else:
        return back_propagate(mode, mother, *parent_probabilities(mode, mm, mf)), \
               back_propagate(mode, father, *parent_probabilities(mode, fm, ff))


def observation_probabilities(mode, observation):
    """Calculate the genotype probabilities for the observation

    Args:
        mode: Mode of inheritance
        observation: The observation in question

    Returns: dict of genotype probabilities
    """
    try:
        a = parent_probabilities(mode, *observation.mother_father)
        b = parent_probabilities_n_children(mode, observation.children, *a)
        probabilities = split_parent_probabilities(b)
    except ZeroDivisionError:
        raise NonMendelianPattern(f'The observations phenotype `{observation}` is not valid given a {mode} '
                                  f'mode of inheritance.')
    else:
        return {FEMALE: probabilities[0], MALE: probabilities[1]}[observation.gender]