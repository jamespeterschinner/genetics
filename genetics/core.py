from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product, chain

from exceptions import NonMendelianPattern

__all__ = ['affected_genotype','ALL_MODES', 'AUTOSOMAL_CHROMOSOMES', 'AUTOSOMAL_DOMINANT', 'AUTOSOMAL_MODES', 'AUTOSOMAL_RECESSIVE',
           'constrain_probabilities', 'count', 'FEMALE', 'FEMALE_SEX_LINKED_CHROMOSOMES', 'FEMALES',
           'genotype_possibilities', 'MALE', 'MALE_SEX_LINKED_CHROMOSOMES', 'MALES', 'normalize_probabilities',
           'phenotypes', 'punnet_occurrences', 'punnet_square', 'SEX_LINKED_MODES', 'VALID_INSTRUCTIONS',
           'VALID_OBSERVATIONS', 'X_LINKED_DOMINANT', 'X_LINKED_RECESSIVE', 'Y_LINKED']

FEMALES = {'f', 'F'}
MALES = {'m', 'M'}
VALID_OBSERVATIONS = set.union(FEMALES, MALES)
VALID_INSTRUCTIONS = {'-', '|'}
FEMALE = 'FEMALE'
MALE = 'MALE'
AUTOSOMAL_DOMINANT = 'AUTOSOMAL_DOMINANT'
AUTOSOMAL_RECESSIVE = 'AUTOSOMAL_RECESSIVE'
X_LINKED_DOMINANT = 'X_LINKED_DOMINANT'
X_LINKED_RECESSIVE = 'X_LINKED_RECESSIVE'
Y_LINKED = 'Y_LINKED'
AUTOSOMAL_CHROMOSOMES = {'aa', 'AA', 'aA', 'Aa'}
FEMALE_SEX_LINKED_CHROMOSOMES = {'xx', 'XX', 'xX', 'Xx'}
MALE_SEX_LINKED_CHROMOSOMES = {'xy', 'XY', 'xY', 'Xy'}

ALL_MODES = (AUTOSOMAL_DOMINANT, AUTOSOMAL_RECESSIVE,
             X_LINKED_DOMINANT, X_LINKED_RECESSIVE,
             Y_LINKED)

AUTOSOMAL_MODES = {AUTOSOMAL_DOMINANT, AUTOSOMAL_RECESSIVE}

SEX_LINKED_MODES = {X_LINKED_DOMINANT, X_LINKED_RECESSIVE,
                    Y_LINKED}


def count(sequence):
    denominator = len(sequence)
    return {key: Fraction(count, denominator) for key, count in Counter(sequence).items()}


def _punnet_square(mother_genotypes, father_genotypes):
    """Combine two possible parent genotypes.
    Equivalent to a punnet square
    """
    return count(list(chain.from_iterable(list(''.join(genotype) for genotype in product(mother, father))
                                          for mother, father in
                                          list(product(mother_genotypes, father_genotypes)))))


def _probabilistic_punnet_square(mother_genotypes, father_genotypes):
    """Combine possible parent genotypes with parent combination likely hood
    """
    result = defaultdict(lambda: 0)
    for (m_genotype, m_probability), (f_genotype, f_probability) in \
            product(mother_genotypes.items(), father_genotypes.items()):
        combination_probability = m_probability * f_probability
        for genotype, probability in _punnet_square(m_genotype, f_genotype).items():
            result[genotype] += combination_probability * probability
    return dict(result)


def punnet_square(mother_genotypes, father_genotypes):
    if isinstance(mother_genotypes, str):
        mother_genotypes = [mother_genotypes]
    if isinstance(father_genotypes, str):
        father_genotypes = [father_genotypes]
    if isinstance(mother_genotypes, list):
        mother_genotypes = count(mother_genotypes)
    if isinstance(father_genotypes, list):
        father_genotypes = count(father_genotypes)

    return _probabilistic_punnet_square(mother_genotypes, father_genotypes)




def affected_genotype(mode, genotype):
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

    elif mode == Y_LINKED and 'Y' in genotype:
        return True
    return False

def phenotypes(genotype, mode):
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

    if observation is not None:
        affected = observation.affected

    if mode == Y_LINKED and gender == FEMALE and affected:
        raise NonMendelianPattern('A female can not be affected if the mode is Y_LINKED')

    if mode in AUTOSOMAL_MODES:
        chromosomes = AUTOSOMAL_CHROMOSOMES
    elif gender:
        chromosomes = {MALE: MALE_SEX_LINKED_CHROMOSOMES,
                       FEMALE: FEMALE_SEX_LINKED_CHROMOSOMES}[gender]
    else:
        chromosomes = MALE_SEX_LINKED_CHROMOSOMES.union(FEMALE_SEX_LINKED_CHROMOSOMES)

    if gender is not None and affected is not None:
        chromosomes = {genotype for genotype in chromosomes if
                       affected_genotype(mode, genotype) == affected}

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
                   if observation in phenotypes(genotype, mode)}
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
               if observation in phenotypes(genotype, mode))
