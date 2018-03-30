# genetics

An educational attempt at modeling mendelian modes of inheritance.

written in pure python.

Pedigrees are input via ASCII art. Eg.
```
f-M
|
m f-m f-M
  |   |
  f   f m M f
```
f = unaffected female;
F = affected female;
m = unaffected male;
M = affected male;
\- = partner;
| = child of

The project aims to generate answers to statistical questions such as:
  what are the first generation possible genotypes given a mode of inheritence?
  
```
>>> p = Pedigree.from_file(r'C:\Users\James\PycharmProjects\genetics\family_tree.txt')
>>> p
                            
  0| |1| |2| |3| |4| |5| |6|
0|f - M                    
 ||                        
1|m   f - m   f - M        
 |    |       |            
2|    f       f   m   M   f
>>> from pprint import pprint
>>> answer = genotypes_given_children(X_LINKED_RECESSIVE, p[0,0])
>>> pprint(answer)
{('XX', 'xY'): Fraction(2, 5),
 ('XX', 'xy'): Fraction(2, 5),
 ('Xx', 'xY'): Fraction(1, 20),
 ('Xx', 'xy'): Fraction(1, 20),
 ('xX', 'xY'): Fraction(1, 20),
 ('xX', 'xy'): Fraction(1, 20)}
 ```
  
