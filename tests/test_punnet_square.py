from fractions import Fraction

from analysis import punnet_square


def test_punnet_square():
    assert punnet_square(['xx', 'xX'], ['xy', 'xY']) == \
           {'xx': Fraction(3, 8), 'xy': Fraction(3, 16), 'xY': Fraction(3, 16), 'Xx': Fraction(1, 8),
            'Xy': Fraction(1, 16), 'XY': Fraction(1, 16)}
