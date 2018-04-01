# genetics

An educational attempt at modeling mendelian modes of inheritance.

written in pure python.

Pedigrees are input via ASCII art. Eg.
```
m-f
  |
  m f-m m   f-m
    |       |
    m m m f M
```
---
```
f = unaffected female
F = affected female
m = unaffected male
M = affected male
- = partner
| = child of
```

The project aims to generate answers to statistical questions such as:
  what are the first generation possible genotypes given a mode of inheritence?
  
```
>>> p = Pedigree.from_file(r'C:\Users\James\PycharmProjects\genetics\family_tree.txt')
>>> p
                            
  0| |1| |2| |3| |4| |5| |6| |7|
0|m - f                        
 |    |                        
1|    m   f - m   m       f - m
 |        |               |    
2|        m   m   m   f   M    

>>> observation_probabilities(X_LINKED_RECESSIVE, p[2,5])
{'xX': Fraction(1, 18), 'XX': Fraction(17, 18)}
>>> observation_probabilities(X_LINKED_RECESSIVE, p[1,6])
{'xX': Fraction(1, 1)}
>>> observation_probabilities(X_LINKED_RECESSIVE, p[1,2])
{'xX': Fraction(1, 9), 'XX': Fraction(8, 9)}
 ```
  
