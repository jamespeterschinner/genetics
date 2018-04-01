import pytest

from exceptions import InvalidObservation
from observation import Observation


def test_observation_raises_error_with_invalid_input():
    with pytest.raises(InvalidObservation):
        Observation('INVALID')


def test_observation_does_not_raise_error_with_valid_input():
    assert Observation('f')
    assert Observation('F')
    assert Observation('m')
    assert Observation('M')


def test_observation_adds_child_to_mother_when_mother_added_to_child():
    a = Observation('f')
    b = Observation('m')
    b.mother = a
    assert b.mother == a
    assert b in a.children

def test_observation_adds_mother_to_child_when_child_added_to_mother():
    a = Observation('f')
    b = Observation('m')
    a.add_children(b)
    assert b in a.children
    assert b.mother == a

def test_observation_adds_father_to_child_when_child_added_to_father():
    a = Observation('m')
    b = Observation('f')
    a.add_children(b)
    assert b in a.children
    assert b.father == a
