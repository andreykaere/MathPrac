#!/usr/bin/python3

from sympy import *
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint

from fractions import Fraction
import numpy as np


def is_point_singular(cubic, point):
    n, x, y, z = symbols('n x y z')
    t = symbols('t')
    
    fx = diff(cubic, x)
    fy = diff(cubic, y)
    fz = diff(cubic, z)

    (x0, y0, z0) = point

    return fx.subs({x: x0, y: y0, z: z0}) == 0 and \
           fy.subs({x: x0, y: y0, z: z0}) == 0 and \
           fz.subs({x: x0, y: y0, z: z0}) == 0 

def get_hessian(c):
    n, x, y, z = symbols('n x y z')
    return Matrix([
       [diff(c, x, x), diff(c, x, y), diff(c, x, z)],
       [diff(c, y, x), diff(c, y, y), diff(c, y, z)],
       [diff(c, z, x), diff(c, z, y), diff(c, z, z)],
    ]).det()


def resultant(f, g, var):
    n, x, y, z = symbols('n x y z')
    matrix = []

    if f == 0 or g == 0:
        return Integer(0)

    f = Poly(f, var)
    g = Poly(g, var)
    n = f.degree()
    m = g.degree()
    size = m + n

    coeff_f = f.all_coeffs()
    coeff_g = g.all_coeffs()

    for i in range(m):
        matrix += [[0]*i + coeff_f + [0]*(size - i - (n+1))]
            
    for i in range(n):
        matrix += [[0]*i + coeff_g + [0]*(size - i - (m+1))]

    return Matrix(matrix).det()


def from_solution_to_rational_points(cubic1, cubic2, solution):
    n, x, y, z, t = symbols('n x y z t')
    xp = solution[0]
    yp = solution[1]
   
    cubic1 = cubic1.subs({x: xp, y: yp}).simplify().expand()

    if cubic1 == 0:
        if cubic2.subs({x: xp, y: yp, z: 1}) == 0:
            return (True, [(xp, yp, 1)])
        else:
            return (True, [])

    if cubic1.as_expr().is_constant():
        return (True, [])
    
    rational_z = get_rational_roots(Poly(cubic1).subs(z, t))

    points = []
    for zp in rational_z:
        if (xp, yp, zp) != (0, 0, 0) and cubic2.subs({x: xp, y: yp, z: zp}) == 0:
            points += [(xp, yp, zp)]

    if points != []: 
        return (True, points)

    return (False, [])    
     
# Considering, that poly is polynomial of variable `t`
def get_rational_roots(poly_t):
    t = symbols('t')
    solutions = set()

    if poly_t.subs(t, 0) == 0:
        solutions.add(0)

    # Getting rid of multiple zero roots so that we can start our enumeration 
    # algorithm 
    while poly_t.subs(t, 0) == 0:
        poly_t = expand(poly_t / t)
    
    if poly_t.as_expr().is_constant():
        return list(solutions)


    coeff_t = Poly(poly_t).all_coeffs()
    last  = coeff_t[-1]
    first = coeff_t[0]
    

    # Function `divisors` gives us positive integer solutions, so we have to
    # consider negative possible ones ourselves
    for j in divisors(last):
        for i in divisors(first):
            if poly_t.subs(t, Fraction(j, i)) == 0:
                solutions.add(Fraction(j, i))
            
            if poly_t.subs(t, Fraction(-j, i)) == 0:
                solutions.add(Fraction(-j, i))

    return list(solutions)
   

def intersection_points(cubic1, cubic2):
    n, x, y, z = symbols('n x y z')
    trans = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]

    a0 = cubic1.coeff(z**3)
    b0 = cubic2.coeff(z**3)

    t = symbols('t')
    res = resultant(cubic1, cubic2, z)

    if res == 0:
        print("Resultant is zero, cubic is reducible, can't proceed ...")
        return (False, [])

    degree = Poly(res).total_degree()

    # Creating set and not array, because we don't care if roots are multiple 
    # or not and in fact don't want to have multiple roots
    solutions = {(0, 0)}

    if res.subs({x: 1, y: 0}) == 0:
        solutions.add((1, 0))


    res_t = Poly(collect(expand(res/y**degree).subs(x/y, t), t), t)
    
    # Reducing our polynom, so that our enumeration algorithm will not take 
    # forever to finish
    coeff_t = res_t.all_coeffs()

    gcd = np.gcd.reduce(coeff_t)
    res_t = simplify(res_t / gcd).expand()

    if res_t.subs(t, 0) == 0:
        solutions.add((0, 1))
        res_t = expand(res_t / t)

    # If there are some other roots, besides 0
    if not res_t.as_expr().is_constant():
        rational_sols = get_rational_roots(Poly(res_t))

        for sol in rational_sols:
            solutions.add((sol.numerator, sol.denominator))
    
   
    points_all = []

    for solution in list(solutions):
        (result, points) = from_solution_to_rational_points(cubic1, cubic2, solution)
        
        if result:
            for point in points:
                points_all += [tuple(Matrix(trans).inv() * Matrix(list(point)))]

    return (True, points_all)
   

def find_inflection_points(cubic):
    n, x, y, z = symbols('n x y z')
    
    hessian = get_hessian(cubic)
    if hessian == 0:
        print("Fatal: Hessian is zero")
        return (False, [])
    
    (res, points) = intersection_points(cubic, hessian)
    
    if not res:
        return (False, [])
    
    if points == []:
        print("Fatal: No rational inflection point exists on the cubic, can't proceed, aborting ...")
        return (False, [])

    return (True, points)


def find_non_singular_inflection_point(cubic):
    (res, points) = find_inflection_points(cubic)

    if not res:
        return (False, ())

    for point in points:
        if not is_point_singular(cubic, point):
            return (True, point)
    
    print("Fatal: There is no rational non-singular inflection point, can't proceed, aborting ...")
    
    return (False, ())


def main():
    print("""Enter your cubic's equation in homogenious coordinates x, y, z:
    For example: x^3 + y^3 + z^3 + 3 x y z
    """)

    n, x, y, z = symbols('n x y z')

    cubic = input()
    cubic = mathematica(cubic)

    (res, point) = find_non_singular_inflection_point(cubic)

    if res:
        (x, y, z) = point
        print("Found rational inflection point: ({} : {} : {})".format(x, y, z))
    else:
        print("There is no rational infleciton point on the cubic")





if __name__ == "__main__":
    main()
