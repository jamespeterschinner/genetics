from core import *

def test_aa_AUTOSOMAL_DOMINANT():
    assert not affected_genotype('aa', AUTOSOMAL_DOMINANT)

def test_aA_AUTOSOMAL_DOMINANT():
    assert affected_genotype('aA', AUTOSOMAL_DOMINANT)

def test_Aa_AUTOSOMAL_DOMINANT():
    assert affected_genotype('Aa', AUTOSOMAL_DOMINANT)

def test_AA_AUTOSOMAL_DOMINANT():
    assert affected_genotype('AA', AUTOSOMAL_DOMINANT)

def test_aa_AUTOSOMAL_RECESSIVE():
    assert affected_genotype('aa', AUTOSOMAL_RECESSIVE)

def test_aA_AUTOSOMAL_RECESSIVE():
    assert not affected_genotype('aA', AUTOSOMAL_RECESSIVE)

def test_Aa_AUTOSOMAL_RECESSIVE():
    assert not affected_genotype('Aa', AUTOSOMAL_RECESSIVE)

def test_AA_AUTOSOMAL_RECESSIVE():
    assert not affected_genotype('AA', AUTOSOMAL_RECESSIVE)

def test_xx_X_LINKED_DOMINANT():
    assert not affected_genotype('xx', X_LINKED_DOMINANT)

def test_xX_X_LINKED_DOMINANT():
    assert affected_genotype('xX', X_LINKED_DOMINANT)

def test_Xx_X_LINKED_DOMINANT():
    assert affected_genotype('Xx', X_LINKED_DOMINANT)

def test_XX_X_LINKED_DOMINANT():
    assert affected_genotype('XX', X_LINKED_DOMINANT)


def test_xy_X_LINKED_DOMINANT():
    assert not affected_genotype('xy', X_LINKED_DOMINANT)


def test_xY_X_LINKED_DOMINANT():
    assert not affected_genotype('xY', X_LINKED_DOMINANT)


def test_Xy_X_LINKED_DOMINANT():
    assert affected_genotype('Xy', X_LINKED_DOMINANT)


def test_XY_X_LINKED_DOMINANT():
    assert affected_genotype('XY', X_LINKED_DOMINANT)


############################

def test_xx_X_LINKED_RECESSIVE():
    assert affected_genotype('xx', X_LINKED_RECESSIVE)


def test_xX_X_LINKED_RECESSIVE():
    assert not affected_genotype('xX', X_LINKED_RECESSIVE)

def test_Xx_X_LINKED_RECESSIVE():
    assert not affected_genotype('Xx', X_LINKED_RECESSIVE)

def test_XX_X_LINKED_RECESSIVE():
    assert not affected_genotype('XX', X_LINKED_RECESSIVE)

def test_xy_X_LINKED_RECESSIVE():
    assert affected_genotype('xy', X_LINKED_RECESSIVE)

def test_xY_X_LINKED_RECESSIVE():
    assert affected_genotype('xY', X_LINKED_RECESSIVE)

def test_Xy_X_LINKED_RECESSIVE():
    assert not affected_genotype('Xy', X_LINKED_RECESSIVE)

def test_XY_X_LINKED_RECESSIVE():
    assert not affected_genotype('XY', X_LINKED_RECESSIVE)

def test_xx_Y_LINKED():
    assert not affected_genotype('xx', Y_LINKED)

def test_xX_Y_LINKED():
    assert not affected_genotype('xX', Y_LINKED)

def test_Xx_Y_LINKED():
    assert not affected_genotype('Xx', Y_LINKED)

def test_XX_Y_LINKED():
    assert not affected_genotype('XX', Y_LINKED)

def test_xy_Y_LINKED():
    assert not affected_genotype('xy', Y_LINKED)

def test_xY_Y_LINKED():
    assert affected_genotype('xY', Y_LINKED)

def test_Xy_Y_LINKED():
    assert not affected_genotype('Xy', Y_LINKED)

def test_XY_Y_LINKED():
    assert affected_genotype('XY', Y_LINKED)

