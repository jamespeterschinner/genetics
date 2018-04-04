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

The project aims to generate answers to statistical questions. Current functionality
is limited. See below for example.
  
```
>>> p2 = Pedigree.from_file(r'C:\Users\James\PycharmProjects\genetics\pedigree_2.txt')
>>> p2
                                  
  0| |1| |2| |3| |4| | | | |5| |6|
0|f - m                          
 ||                              
1|m   f - m   m             f - m
 |    |                     |    
2|    m   m   m   f         M    
>>> genotypes_n_mode(X_LINKED_RECESSIVE, p2[0,0])
{'Xx': Fraction(288, 823), 'xX': Fraction(535, 823)}
>>> genotypes_n_mode(X_LINKED_RECESSIVE, p2[1,1])
{'XX': Fraction(8, 9), 'xX': Fraction(1, 9)}
>>> genotypes_n_mode(X_LINKED_RECESSIVE, p2[2,4])
{'XX': Fraction(17, 18), 'xX': Fraction(1, 18)}
 ```
  
