from analysis import phenotype_given_genotype
from constants import *


def test_aa_AUTOSOMAL_DOMINANT():
    assert phenotype_given_genotype('aa', AUTOSOMAL_DOMINANT) == 'fm'


def test_aA_AUTOSOMAL_DOMINANT():
    assert phenotype_given_genotype('aA', AUTOSOMAL_DOMINANT) == 'FM'


def test_Aa_AUTOSOMAL_DOMINANT():
    assert phenotype_given_genotype('Aa', AUTOSOMAL_DOMINANT) == 'FM'


def test_AA_AUTOSOMAL_DOMINANT():
    assert phenotype_given_genotype('AA', AUTOSOMAL_DOMINANT) == 'FM'


def test_aa_AUTOSOMAL_RECESSIVE():
    assert phenotype_given_genotype('aa', AUTOSOMAL_RECESSIVE) == 'FM'


def test_aA_AUTOSOMAL_RECESSIVE():
    assert phenotype_given_genotype('aA', AUTOSOMAL_RECESSIVE) == 'fm'


def test_Aa_AUTOSOMAL_RECESSIVE():
    assert phenotype_given_genotype('Aa', AUTOSOMAL_RECESSIVE) == 'fm'


def test_AA_AUTOSOMAL_RECESSIVE():
    assert phenotype_given_genotype('AA', AUTOSOMAL_RECESSIVE) == 'fm'


def test_xx_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('xx', X_LINKED_DOMINANT) == 'f'


def test_xX_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('xX', X_LINKED_DOMINANT) == 'F'


def test_Xx_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('Xx', X_LINKED_DOMINANT)  == 'F'


def test_XX_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('XX', X_LINKED_DOMINANT) == 'F'


def test_xy_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('xy', X_LINKED_DOMINANT) == 'm'


def test_xY_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('xY', X_LINKED_DOMINANT) == 'm'


def test_Xy_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('Xy', X_LINKED_DOMINANT) == 'M'


def test_XY_X_LINKED_DOMINANT():
    assert phenotype_given_genotype('XY', X_LINKED_DOMINANT) == 'M'


############################

def test_xx_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('xx', X_LINKED_RECESSIVE) == 'F'


def test_xX_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('xX', X_LINKED_RECESSIVE) == 'f'


def test_Xx_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('Xx', X_LINKED_RECESSIVE) == 'f'


def test_XX_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('XX', X_LINKED_RECESSIVE) == 'f'


def test_xy_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('xy', X_LINKED_RECESSIVE) == 'M'


def test_xY_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('xY', X_LINKED_RECESSIVE) == 'M'


def test_Xy_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('Xy', X_LINKED_RECESSIVE) == 'm'


def test_XY_X_LINKED_RECESSIVE():
    assert phenotype_given_genotype('XY', X_LINKED_RECESSIVE) == 'm'


def test_xx_Y_LINKED():
    assert phenotype_given_genotype('xx', Y_LINKED) == 'f'


def test_xX_Y_LINKED():
    assert phenotype_given_genotype('xX', Y_LINKED) == 'f'


def test_Xx_Y_LINKED():
    assert phenotype_given_genotype('Xx', Y_LINKED) == 'f'


def test_XX_Y_LINKED():
    assert phenotype_given_genotype('XX', Y_LINKED) == 'f'


def test_xy_Y_LINKED():
    assert phenotype_given_genotype('xy', Y_LINKED) == 'M'


def test_xY_Y_LINKED():
    assert phenotype_given_genotype('xY', Y_LINKED) == 'M'


def test_Xy_Y_LINKED():
    assert phenotype_given_genotype('Xy', Y_LINKED) == 'M'


def test_XY_Y_LINKED():
    assert phenotype_given_genotype('XY', Y_LINKED) == 'M'
