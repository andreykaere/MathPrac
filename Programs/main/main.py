#!/usr/bin/python3

from rational_points.point_sum import (FindRational, RevNFind)
from rational_points.draw_sum import DrawSumN
from weierstrass_form.weierstrass_main import weierstrass_form

from fractions import Fraction
from sympy import *
from sympy.parsing.mathematica import mathematica
import numpy as np

n, x, y, z = symbols('n x y z')
cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"
cubic = mathematica(cubic)
cubic = cubic.subs(n, 4)

(res, ((a, b), trans)) = weierstrass_form(cubic)

inv_trans = np.array(list(map(lambda i: Fraction(i), Matrix(trans).inv())))
inv_trans = np.reshape(inv_trans, (3, 3))

R = 700
P = FindRational(a, b, R)
if P[1] == 1:
    print("Start point is :", P[0])
    x = P[0][0]
    y = P[0][1]
   
    n = 50
    RevNFind(a, b, x, y, n, inv_trans)
    
    n = 4
    DrawSumN(a, b, x, y, n)    
    
else:
    print("No good rational points found, can't apply algoritm")

