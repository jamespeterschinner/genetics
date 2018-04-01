from analysis import normalize_probabilities
from fractions import Fraction

def test_normalize_probabilities():
    assert normalize_probabilities({1: Fraction(1, 4), 2: Fraction(2, 3)}) == \
           {1: Fraction(3, 11), 2: Fraction(8, 11)}
