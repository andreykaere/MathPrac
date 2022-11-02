#!/usr/bin/python3

from sympy import *
# from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.mathematica import mathematica

from fractions import Fraction
import numpy as np
import math

# We assume, that cubic is given in coordinates x = x1, y = x2, z = x3
# n, x, y, z = symbols('n x y z')
# n = symbols('n', integer = True)

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



# Map point `point` to (0 : 1 : 0)
def weirstrass_form_step1(cubic, point):
    n, x, y, z = symbols('n x y z')
    (xp, yp, zp) = point

    matrix0 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]

    if zp != 0:
        if xp == 0:
            (zp, xp)  = (xp, zp)
            matrix0 = [
                [0, 0, 1],
                [0, 1, 0],
                [1, 0, 0],
            ]

        if yp == 0:
            (zp, yp)  = (yp, zp)
            matrix0 = [
                [1, 0, 0],
                [0, 0, 1],
                [0, 1, 0],
            ]
    
    matrix1 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    
    if xp == 0:
        (xp, yp)  = (yp, xp)
        matrix1 = [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
        ]

    matrix2 = [
        [Fraction(yp, xp), -1, 1],
        [Fraction(1, xp), 0, 0],
        [0, 0, 1],
    ]

    matrix = np.dot(matrix2,
             np.dot(matrix1, 
                    matrix0
             ).tolist()
             ).tolist()
    
    a, b, c = symbols('a b c')

    (x1, y1, z1) = tuple(Matrix(matrix).inv() * Matrix([a, b, c]))

    cubic = simplify(cubic.subs({x: x1, y: y1, z: z1}))
    cubic = cubic.subs({a: x, b: y, c: z})
    
    # print(cubic)

    return (matrix, cubic)

    

def weirstrass_form_step2(cubic):
    n, x, y, z = symbols('n x y z')
    c11, c31, c21, c23, c13, c33 = symbols('c11 c31 c21 c23 c13 c33') 

    matrix = Matrix([
        [c11, 0, c13],
        [c21, 1, c23],
        [c31, 0, c33],
    ])

    p = diff(cubic, x).subs({x: 0, y: 1, z: 0})
    q = diff(cubic, z).subs({x: 0, y: 1, z: 0})

    # Matrix([p, q]) * (Matrix([
    #     [c11, c13],
    #     [c21, c23],
    #     [c31, c33],
    # ]).T)

    print(p)
    print(q)

    if p == 0 and q == 0:
        print("SOMETHING IS WRONG: p == 0 and q == 0")

    

    if p != 0:
        matrix.subs({
            c11: 1,
            c13: -Fraction(q, p),
            c21: 1,
            c23: -Fraction(q, p),
            c31: 1,
            c33: -Fraction(q, p) - 1,
        })
    
    if q != 0:
        matrix.subs({
            c13: 1,
            c11: -Fraction(p, q),
            c23: 1,
            c21: -Fraction(p, q),
            c33: 1,
            c31: -Fraction(p, q) - 1,
        })

    # solve([
    # (p * c11 + q * c13) == 0,
    # (p * c21 + q * c23) == 0,
    # (p * c31 + q * c33) != 0,
    # ], [c11, c13, c21, c23,])
    # solve_rational_inequalities([[
    #     ((Poly(p * c11 + q * c31), Poly(0, c11)), '=='),
    #     ((Poly(p * c21 + q * c23), Poly(0, c21)), '=='),
    #     ((Poly(p * c31 + q * c33), Poly(0, c31)), '!='),
    # ]])

    print(matrix)



def weirstrass_form_step3(cubic):
    n, x, y, z = symbols('n x y z')
    pass
    
def weirstrass_form(cubic):
    n, x, y, z = symbols('n x y z')

    (points, trans0) = points_of_inflection(cubic)
    point = points[0]
    
    (trans1, cubic) = weirstrass_form_step1(cubic, point)
    (trans2, cubic) = weirstrass_form_step2(cubic)
    (trans3, cubic) = weirstrass_form_step3(cubic)

   
    trans = np.dot(trans3,
            np.dot(trans2, 
            np.dot(trans1,
                   trans0
            ).tolist()
            ).tolist()
            ).tolist() 

    return ((a, b), trans)

    





def main():
    # print("Enter your cubic's equation in homogenious coordinates x, y, z:")

    # cubic = input()
    n, x, y, z = symbols('n x y z')

    cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"

    cubic = mathematica(cubic)

    cubic = cubic.subs(n, 4)

    # weirstrass_form(cubic)

    print(weirstrass_form_step1(cubic, (-1, 1, 0)))
    (trans, cubic) = weirstrass_form_step1(cubic, (-1, 1, 0))
    print(weirstrass_form_step2(cubic))



main()

# print(resultant(get_gessian(cubic), cubic, z))
