from constants import *
from pedigree import Pedigree
from analysis import *
from observation import Observation
from pprint import pprint as pp

p1 = Pedigree.from_file(r'C:\Users\James\PycharmProjects\genetics\pedigree_1.txt')
p2 = Pedigree.from_file(r'C:\Users\James\PycharmProjects\genetics\pedigree_2.txt')
mode = X_LINKED_RECESSIVE