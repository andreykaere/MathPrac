#!/usr/bin/python3

from sympy import *
# from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint
# from sympy.ntheory import primefactors

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

    # # Just a test 
    # coeff_f = [i for i in range(1, n+2)]
    # coeff_g = [-i for i in range(1, m+2)]

    for i in range(m):
        matrix += [[0]*i + coeff_f + [0]*(size - i - (n+1))]
            
    for i in range(n):
        matrix += [[0]*i + coeff_g + [0]*(size - i - (m+1))]

    return Matrix(matrix).det()


def from_solution_to_rational_points(cubic1, cubic2, solution):
    n, x, y, z, t = symbols('n x y z t')
    # x_, y_, z_ = symbols('n x y z')
    xp = solution[0]
    yp = solution[1]
   
    # print("point", (xp, yp))
    # print("cubic1 before point", cubic1)
    cubic1 = cubic1.subs({x: xp, y: yp}).simplify().expand()
    # print("cubic1 after point", cubic1)

    if cubic1 == 0:
        return (True, [(xp, yp, 1)])

    if cubic1.as_expr().is_constant():
        return (True, [])
    
    # print("foo")
    # print("rational_z, cubic1:", cubic1)
    # print("rational_z, cubic1:", Poly(cubic1).subs(z, t))
    # print("bar")
    rational_z = get_rational_roots(Poly(cubic1).subs(z, t))

    points = []
    for zp in rational_z:
        if (xp, yp, zp) != (0, 0, 0) and cubic2.subs({x: xp, y: yp, z: zp}) == 0:
            points = [(xp, yp, zp)]

    if points != []: 
        return (True, points)

    return (False, [])    
     
    
def fix_zero_leading_coefficients(cubic, hessian):
    n, x, y, z = symbols('n x y z')
    
    flag = 0
    # Searching for point that does not lie on both cubics
    for i in range(-2, 2):
        for j in range(-2, 2):
            if i != 0 and j != 0 and cubic.subs({x: i, y: j, z: 1}) != 0 and \
               hessian.subs({x: i, y: j, z: 1}) != 0:
                print("I found point, that does not lie on both!")
                flag = 1
                break
        
        if flag == 1:
            break

    if flag == 0:
        print("Cubic and hessian coincide!")
   
    point = (i, j, 1)
    
    matrix = [
        [1, 0, i],
        [0, 1, j],
        [0, 0, 1],
    ]
    matrix = Matrix(matrix).inv().tolist()
    
    a, b, c = symbols('a b c')
    (x1, y1, z1) = tuple(Matrix(matrix).inv() * Matrix([a, b, c]))
    
    cubic = simplify(cubic.subs({x: x1, y: y1, z: z1})).expand()
    cubic = cubic.subs({a: x, b: y, c: z})
    hessian = simplify(hessian.subs({x: x1, y: y1, z: z1})).expand()
    hessian = hessian.subs({a: x, b: y, c: z})
    
    return (cubic, hessian, matrix)



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
        # print("poly_t is constant")
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

    # TODO: Not sure if we need this
    # if a0 * b0 == 0:
    #     (cubic1, cubic2, trans) = fix_zero_leading_coefficients(cubic1, cubic2)

    # print("cubic1, after fix", cubic1)
    # print("cubic2, after fix", cubic2)

    t = symbols('t')
    res = resultant(cubic1, cubic2, z)

    # print("This is resultant", res)

    if res == 0:
        print("Resultant is zero, cubic is reducible, can't proceed ...")
        return (False, [])


    degree = Poly(res).total_degree()

    if degree < 9:
        print("WARNING: Resultant is degenerated")
        # return 
    
    # Creating set and not array, because we don't care if roots are multiple 
    # or not and in fact don't want to have multiple roots
    solutions = {(0, 0)}

    if res.subs({x: 1, y: 0}) == 0:
        solutions.add((1, 0))


    res_t = Poly(collect(expand(res/y**degree).subs(x/y, t), t), t)
    
    # Reducing our polynom, so that our enumeration algorithm will not take 
    # forever to finish
    coeff_t = res_t.all_coeffs()

    # print(coeff_t)
    # print("this is res_t", res_t)
    gcd = np.gcd.reduce(coeff_t)
    res_t = simplify(res_t / gcd).expand()


    if res_t.subs(t, 0) == 0:
        solutions.add((0, 1))
        res_t = expand(res_t / t)
    
    rational_sols = get_rational_roots(Poly(res_t))

    # print(rational_sols)

    for sol in rational_sols:
        solutions.add((sol.numerator, sol.denominator))


    # print(solutions)
    
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
    n, x, y, z = symbols('n x y z')

    # cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"
    # cubic = "y^2 z - x^3 - x^2 z"
    # cubic = "-x^3 - x^2*z + y^2*z + 2*y*z^2 + z^3"
    # cubic = "-x^3 - 4*x^2*z + y^2*z - 5*x*z^2 - 2*z^3"
    # cubic = "-x^3 - 4*x^2*z + y^2*z - 5*x*z^2 + 2*y*z^2 - z^3"
    # cubic = "-x^3 - 3*x^2*z + y^2*z - 3*x*z^2"

    # cubic = "5 y^3 + z^2 x + y^2 x - 34 y^2 z"
    # cubic = "(x - y) (y^2 - x^2 + z x) - x^2 y"
    cubic = "(x - y) (y^2 - x^2 + z x) - x^2 y"
    cubic = "-x^3 + x*y^2 - y^3 + x^2*z - x*y*z"
    # cubic = "(x - z) (x z - y^2)"

    # cubic = "x^3*z + x*y^2*z + x^2*z^2 + y^2*z^2"
    cubic = mathematica(cubic)
    # cubic = cubic.subs(n, 4)

    print("Hessian", get_hessian(cubic))
    # print(cubic)
    # print(get_hessian(cubic))

    print(find_non_singular_inflection_point(cubic))
    # print(cubic)



if __name__ == "__main__":
    main()
