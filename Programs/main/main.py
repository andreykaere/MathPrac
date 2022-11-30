#!/usr/bin/python3

from RationalPoints.PointSum import *
from RationalPoints.DrawSum import *
from WeierstrassForm.WeierstrassForm import weierstrass_form

from sympy import *
from sympy.parsing.mathematica import mathematica
import numpy as np

n, x, y, z = symbols('n x y z')
cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"
cubic = mathematica(cubic)
cubic = cubic.subs(n, 6)
# cubic = cubic.subs(n, 8)

(res, ((a, b), trans)) = weierstrass_form(cubic)


print(Matrix(trans) * Matrix([-1, 0, 1]))
print(Matrix(trans) * Matrix([-1, 1, 0]))
print(Matrix(trans) * Matrix([0, 1, -1]))
# print(Matrix(trans).inv() * Matrix([246, 2106, 1]))

# print(Matrix(trans))
# print(Matrix(trans).inv())
inv_trans = np.array(list(map(lambda i: Fraction(i), Matrix(trans).inv())))
inv_trans = np.reshape(inv_trans, (3, 3))

# print(inv_trans)
# a = -3260115
# b = 2265609582
# print(a, b)
# print()
#-----------------------
R = 700

P = FindRational(a, b, R)
if P[1] == 1:
    print("Start point is :", P[0])
    x = P[0][0]
    y = P[0][1]
   
    # # Special point for a picture for N = 4
    # x = Fraction(246, 1)
    # y = Fraction(2106, 1)     
    
    n = 50
    
    RevNFind(a, b, x, y, n, inv_trans)
    
    n = 4
    DrawSumN(a, b, x, y, n)    

    
else:
    print("No good rational points found, can't apply algoritm")





# x = Fraction(-573, 1)
# y = Fraction(7020, 1)

#x = Fraction(246, 1)
#y = Fraction(2106, 1)

#-----------------------
