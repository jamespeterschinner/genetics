from fractions import Fraction
import pytest
from exceptions import InvalidState
from analysis import genotype_possibilities
from constants import *


def test_genotype_possibilities():
    assert genotype_possibilities(AUTOSOMAL_DOMINANT, MALE) == \
           {'Aa': Fraction(1, 4), 'AA': Fraction(1, 4), 'aA': Fraction(1, 4), 'aa': Fraction(1, 4)}

    assert genotype_possibilities(AUTOSOMAL_DOMINANT, FEMALE) == \
           {'Aa': Fraction(1, 4), 'AA': Fraction(1, 4), 'aA': Fraction(1, 4), 'aa': Fraction(1, 4)}

    assert genotype_possibilities(AUTOSOMAL_DOMINANT, FEMALE, True) == \
           {'Aa': Fraction(1, 3), 'AA': Fraction(1, 3), 'aA': Fraction(1, 3)}

    assert genotype_possibilities(AUTOSOMAL_DOMINANT, MALE, True) == \
           {'Aa': Fraction(1, 3), 'AA': Fraction(1, 3), 'aA': Fraction(1, 3)}

    assert genotype_possibilities(AUTOSOMAL_DOMINANT, FEMALE, False) == \
           {'aa': Fraction(1, 1)}

    assert genotype_possibilities(AUTOSOMAL_DOMINANT, MALE, False) == \
           {'aa': Fraction(1, 1)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, FEMALE) == \
           {'Xx': Fraction(1, 4), 'xX': Fraction(1, 4), 'XX': Fraction(1, 4), 'xx': Fraction(1, 4)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, MALE) == \
           {'Xy': Fraction(1, 4), 'xY': Fraction(1, 4), 'xy': Fraction(1, 4), 'XY': Fraction(1, 4)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, FEMALE, True) == \
           {'xx': Fraction(1, 1)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, FEMALE, False) == \
           {'xX': Fraction(1, 3), 'Xx': Fraction(1, 3), 'XX': Fraction(1, 3)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, MALE, True) == \
           {'xy': Fraction(1, 2), 'xY': Fraction(1, 2)}

    assert genotype_possibilities(X_LINKED_RECESSIVE, MALE, False) == \
           {'XY': Fraction(1, 2), 'Xy': Fraction(1, 2)}

    assert genotype_possibilities(X_LINKED_DOMINANT, FEMALE) == \
           {'Xx': Fraction(1, 4), 'xX': Fraction(1, 4), 'XX': Fraction(1, 4), 'xx': Fraction(1, 4)}

    assert genotype_possibilities(X_LINKED_DOMINANT, MALE) == \
           {'Xy': Fraction(1, 4), 'xY': Fraction(1, 4), 'xy': Fraction(1, 4), 'XY': Fraction(1, 4)}

    assert genotype_possibilities(X_LINKED_DOMINANT, FEMALE, True) == \
           {'xX': Fraction(1, 3), 'Xx': Fraction(1, 3), 'XX': Fraction(1, 3)}

    assert genotype_possibilities(X_LINKED_DOMINANT, FEMALE, False) == \
           {'xx': Fraction(1, 1)}

    assert genotype_possibilities(X_LINKED_DOMINANT, MALE, True) == \
           {'XY': Fraction(1, 2), 'Xy': Fraction(1, 2)}

    assert genotype_possibilities(X_LINKED_DOMINANT, MALE, False) == \
           {'xy': Fraction(1, 2), 'xY': Fraction(1, 2)}

    assert genotype_possibilities(Y_LINKED, FEMALE) == \
           {'Xx': Fraction(1, 4), 'xX': Fraction(1, 4), 'XX': Fraction(1, 4), 'xx': Fraction(1, 4)}

    assert genotype_possibilities(Y_LINKED, FEMALE) == genotype_possibilities(Y_LINKED, FEMALE, False)

    with pytest.raises(InvalidState):
        genotype_possibilities(Y_LINKED, FEMALE, True)