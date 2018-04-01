from itertools import product, chain
from fractions import Fraction
from collections import Counter, defaultdict

def count(sequence):
    denominator = len(sequence)
    return {key: Fraction(count, denominator) for key, count in Counter(sequence).items()}

def _parent_combinations(mother_genotypes, father_genotypes):
    """Return the parent combinations for two lists of parents
    """
    return list(product(mother_genotypes, father_genotypes))


def _punnet_square(mother_genotypes, father_genotypes):
    """Combine two possible parent genotypes.
    Equivalent to a punnet square
    """
    return count(list(chain.from_iterable(list(''.join(genotype) for genotype in product(mother, father))
                                         for mother, father in
                                         _parent_combinations(mother_genotypes, father_genotypes))))

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
