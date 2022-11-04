#!/usr/bin/python3

from sympy import *
# from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint
# from sympy.ntheory import primefactors

from fractions import Fraction
import numpy as np



def get_gessian(c):
    n, x, y, z = symbols('n x y z')
    return Matrix([
       [diff(c, x, x), diff(c, x, y), diff(c, x, z)],
       [diff(c, y, x), diff(c, y, y), diff(c, y, z)],
       [diff(c, z, x), diff(c, z, y), diff(c, z, z)],
    ]).det()


def resultant(f, g, var):
    n, x, y, z = symbols('n x y z')
    matrix = []

    f = Poly(f, var)
    g = Poly(g, var)
    
    n = f.degree()
    m = g.degree()
    size = m + n

    coeff_f = f.all_coeffs()
    coeff_g = g.all_coeffs()

    # coeff_f = [i for i in range(1, n+2)]
    # coeff_g = [-i for i in range(1, m+2)]

    for i in range(m):
        matrix += [[0]*i + coeff_f + [0]*(size - i - (n+1))]
            
    for i in range(n):
        matrix += [[0]*i + coeff_g + [0]*(size - i - (m+1))]

    return Matrix(matrix).det()


def from_solution_to_point(cubic, sol):
    n, x, y, z = symbols('n x y z')
    if sol == Fraction(0):
        pass
        # solutions_xyz += [(0, 1, )]
        # cubic 
    
    if sol == Fraction(j, i):
        pass


def fix_zero_leading_coefficients(cubic, gessian):
    n, x, y, z = symbols('n x y z')
    pass


def points_of_inflection(cubic):
    n, x, y, z = symbols('n x y z')
    trans = []
    gessian = get_gessian(cubic)

    # if a0 b0 = 0, fix it 
    a0 = Poly(cubic, z).all_coeffs()[0]
    b0 = Poly(gessian, z).all_coeffs()[0]

    if a0 * b0 == 0:
        (cubic, gessian, trans) = fix_zero_leading_coefficients(cubic, gessian)


    t = symbols('t')
    res = resultant(cubic, gessian, z)
    res_t = Poly(collect(expand(res/y**9).subs(x/y, t), t), t)
    
    # Reducing our polynom, so that our enumeration algorithm will not take 
    # forever to finish
    coeff_t = res_t.all_coeffs()
    gcd = np.gcd.reduce(coeff_t)
    res_t = simplify(res_t / gcd)

    # Creating set and not array, because we don't care if roots are multiple 
    # or not
    solutions = set()

    if res_t.subs(t, 0) == 0:
        solutions.add(Fraction(0))

    # Getting rid of multiple zero roots so that we can start our enumeration 
    # algorithm 
    while res_t.subs(t, 0) == 0:
        res_t = expand(res_t / t)
    
    # Updating data after reducing
    coeff_t = res_t.all_coeffs()
    last  = coeff_t[-1]
    first = coeff_t[0]

    for j in divisors(last):
        for i in divisors(first):
            if res_t.subs(t, Fraction(j, i)) == 0:
                solutions.add(Fraction(j, i))
            
            if res_t.subs(t, Fraction(-j, i)) == 0:
                solutions.add(Fraction(-j, i))
    
    points = []

    for solution in list(solutions):
        point = from_solution_to_point(solution)
        points.append(point)

    # TODO: printing when there are no rational points on cubic

    return (points, trans)
