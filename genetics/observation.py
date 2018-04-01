from exceptions import InvalidObservation
from constants import *
from analysis import *

class Observation(str):

    @property
    def gender(self):
        gender = None
        if self in FEMALES:
            gender = FEMALE
        elif self in MALES:
            gender = MALE
        return gender

    ### PARTNER ###
    @property
    def partner(self):
        return self._partner

    @partner.setter
    def partner(self, partner):
        self._partner = partner
        if partner is not None:
            partner._partner = self

    ### MOTHER ###
    @property
    def mother(self):
        return self._mother

    @mother.setter
    def mother(self, mother):
        self._mother = mother
        if mother is not None:
            mother.add_children(self)

    ### FATHER ###
    @property
    def father(self):
        return self._father

    @father.setter
    def father(self, father):
        self._father = father
        if father is not None:
            father.add_children(self)

    ### PARENT ###
    @property
    def parent(self):
        result = None
        if self.mother:
            result = self.mother
        elif self.father:
            result = self.mother
        return result

    ### CHILDREN ###
    @property
    def children(self):
        return self._children

    ### SIBLINGS ###
    @property
    def siblings(self):
        children = ()
        if self.mother:
            children = self.mother.children
        elif self.father:
            children = self.father.children
        return children

    ### ROOT ###
    @property
    def root(self):
        return self.mother is None and self.father is None

    ### AFFECTED ###
    @property
    def affected(self):
        return self.isupper()

    ### MOTHER FATHER ###
    @property
    def mother_father(self):
        return {True: (self, self.partner),
                False: (self.partner, self)}[self.gender == FEMALE]

    def add_children(self, value):
        # Validate children are Observations
        if isinstance(value, (tuple, list)):
            if all(map(lambda x: isinstance(x, Observation), value)):
                self._children += value
        elif isinstance(value, Observation):
            self._children += (value,)

        mother, father = self.mother_father
        for child in self._children:
            child._mother, child._father = mother, father



    def __new__(cls, value):
        if value not in VALID_OBSERVATIONS:
            raise InvalidObservation
        return super().__new__(cls, value)

    def __init__(self, value):
        super().__init__()
        self._mother = None
        self._father = None
        self._partner = None
        self._children = ()
        self._index = None


