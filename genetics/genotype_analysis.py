from collections import defaultdict
from functools import reduce
from itertools import product
from operator import mul

from core import FEMALE, MALE
from core import constrain_probabilities
from core import genotype_possibilities
from core import normalize_probabilities
from core import punnet_occurrences
from core import punnet_square
from exceptions import NonMendelianPattern

__all__ = ['parent_likelihood', 'genotypes_n_mode']

def _forward_propagate(mode, child, genotype_probabilities):
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


def parent_likelihood(mode, children, mother_genotypes, father_genotypes):
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
            parent_likelihood(mode, child.children,
                              *_forward_propagate(mode, child, genotypes_probabilities))
        )
    return reduce(mul, result)


def _fit_parent_probabilities(mode, children, mother_genotypes, father_genotypes):
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
            parent_likelihood(mode, children, m_genotype, f_genotype) * (f_prob * m_prob)
        for (m_genotype, m_prob), (f_genotype, f_prob) in
        product(mother_genotypes.items(), father_genotypes.items())
    })


def _parent_to_children_probabilities(parent_probabilities):
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


def _split_parent_probabilities(probabilities):
    mother = defaultdict(lambda: 0)
    father = defaultdict(lambda: 0)
    for (m_genotype, f_genotype), probability in probabilities.items():
        mother[m_genotype] += probability
        father[f_genotype] += probability
    return dict(mother), dict(father)


def _back_propagate(mode, observation, mother_genotypes, father_genotypes):
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
    a = _fit_parent_probabilities(mode, observation.siblings, mother_genotypes, father_genotypes)
    b = _parent_to_children_probabilities(a)
    return constrain_probabilities(mode, observation, b)


def _parent_probabilities(mode, mother, father):
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
        return mother_genotypes, _back_propagate(mode, father, *_parent_probabilities(mode, fm, ff))
    elif father is None:
        return _back_propagate(mode, mother, *_parent_probabilities(mode, mm, mf)), father_genotypes
    else:
        return _back_propagate(mode, mother, *_parent_probabilities(mode, mm, mf)), \
               _back_propagate(mode, father, *_parent_probabilities(mode, fm, ff))


def genotypes_n_mode(mode, observation):
    """Calculate the genotype probabilities for the observation

    Args:
        mode: Mode of inheritance
        observation: The observation in question

    Returns: dict of genotype probabilities
    """
    try:
        a = _parent_probabilities(mode, *observation.mother_father)
        b = _fit_parent_probabilities(mode, observation.children, *a)
        probabilities = _split_parent_probabilities(b)
    except ZeroDivisionError:
        raise NonMendelianPattern(f'The observations phenotype `{observation}` is not valid given a {mode} '
                                  f'mode of inheritance.')
    else:
        return {FEMALE: probabilities[0], MALE: probabilities[1]}[observation.gender]
