from PointSum3 import *
from DrawSum import *

#-----------------------
N = 4
a = -432 * (N**4) - 2592 * (N**3) - 3240 * (N**2) + 4536 * N + 7533
b = 3456 * (N**6) + 31104 * (N**5) + 85536 * (N**4) + 15552 * (N**3) - 250776 * (N**2) - 239112 * N + 68526
x = Fraction(246, 1)
y = Fraction(2106, 1)


n = 9

RevNFind(a, b, x, y, n, N)

n = 9
DrawSumN(a, b, x, y, n)

#-----------------------