from PointSum3 import *
from DrawSum import *
from WeierstrassForm import *

n, x, y, z = symbols('n x y z')

cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"

cubic = mathematica(cubic)

cubic = cubic.subs(n, 4)

((a, b), trans) = weirstrass_form(cubic)

print(a, b)

#-----------------------
N = 4
N = Fraction(4)


x = Fraction(246, 1)
y = Fraction(2106, 1)


n = 30

RevNFind(a, b, x, y, n, trans)

n = 9
DrawSumN(a, b, x, y, n)

#-----------------------
