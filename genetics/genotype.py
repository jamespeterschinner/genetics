from itertools import product

class Genotype(str):

    def __add__(self, other):
        female, male = sorted([self, other], key=lambda x: 'y' in x.lower())
        return [self.__class__(''.join(i)) for i in product(female, male)]
