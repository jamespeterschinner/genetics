from constants import *
from itertools import product, chain
from fractions import Fraction
from functools import partial, reduce
from operator import mul

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


def valid_genotypes(mode, gender, affected=None, observation=None):
    """Return valid genotypes for a given observation and mode

    Args:
        mode: The mode of inheritance
        gender: Optionally supply the gender
        affected: True the genotypes will be affected for the given mode

    Returns: list of valid genotypes
    """
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

    return chromosomes


def parent_combinations(mother_genotypes, father_genotypes):
    return list(product(mother_genotypes, father_genotypes))

def punnet_square(mother_genotypes, father_genotypes):
    """Combine two possible parent genotypes.
    Equivalent to a punnet square"""
    return list(chain.from_iterable(list(''.join(genotype) for genotype in product(mother, father))
                                    for mother, father in
                                    parent_combinations(mother_genotypes, father_genotypes)))

def constrain_punnet(mode, observation, punnet):
    return [genotype for genotype in punnet if
            observation in phenotype_given_genotype(genotype, mode)]

def punnet_occurrences(mode, observation, punnet):
    """Calculate the times the observation would have been made given a punnet square"""
    try:
        return Fraction(len(constrain_punnet(mode, observation, punnet)), len(punnet))
    except ZeroDivisionError:
        return 0


def next_punnet(mode, observation, punnet):
    """Combine a punnet square with an observation

    The idea is to constrain the supplied punnet with the observation
    and return a new punnet square with the possible partner genotypes

    Args:
        mode: Mode of inheritance
        observation: The observation to constrain the supplied punnet
        punnet: A list of possible genotypes

    Returns: list of possible genotypes
    """
    mother, father = observation.mother_father

    if mother is observation:
        mother_genotypes = constrain_punnet(mode, mother, punnet)
    elif father is observation:
        father_genotypes = constrain_punnet(mode, father, punnet)

    if mother is None:
        mother_genotypes = valid_genotypes(mode, FEMALE)
    elif father is None:
        father_genotypes = valid_genotypes(mode, MALE)
    elif mother is not observation:
        mother_genotypes = valid_genotypes(mode, FEMALE, mother.affected)
    elif father is not observation:
        father_genotypes = valid_genotypes(mode, MALE, father.affected)

    return punnet_square(mother_genotypes, father_genotypes)



def punnet_probability(mode, observation, punnet):
    """Calculate the probability of children observations

    Args:
        mode: Mode of inheritance
        children: An iterable of children observations
        punnet_func: A callable that returns a list of possible genotypes

    Returns: Fraction indicating probability
    """
    probabilities = []
    for child in observation.children:
        probabilities.append(
            punnet_occurrences(mode, child, punnet)
        )
        probabilities.append(
            punnet_probability(mode, child, next_punnet(mode, child, punnet))
        )

    return reduce(mul, probabilities, 1)


def normalize_probabilities(probabilities):
    prior = Fraction(1, len(probabilities))
    normalize = sum([i * prior for i in probabilities.values()])
    return {key: (value * prior) / normalize for key, value in probabilities.items()
            if value != 0}

def genotypes_given_children(mode, observation):
    """Constrain the observations genotype from the observations
    children.

    Args:
        mode: The mode of inheritance
        observation: The observations to constrain.
    """

    mother, father = observation.mother_father
    mother_genotypes = valid_genotypes(mode, FEMALE, observation=mother)
    father_genotypes = valid_genotypes(mode, MALE, observation=father)
    parent_genotype_probabilities = {}
    for mother_genotype, father_genotype in product(mother_genotypes, father_genotypes):
        parent_genotype_probabilities.update(
            {(mother_genotype, father_genotype):
            punnet_probability(mode, observation, punnet_square(mother_genotype, father_genotype))}
        )
    return normalize_probabilities(parent_genotype_probabilities)
