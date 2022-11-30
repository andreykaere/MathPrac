#!/usr/bin/python3

from sympy import *
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint

from fractions import Fraction
import numpy as np

# import importlib
# from weierstrass_form.inflection_points import find_non_singular_inflection_point
# from inflection_points import find_non_singular_inflection_point
# importlib.import_module('inflection_points')
# importlib.import_module('weierstrass_form.inflection_points')
# from weierstrass_form.inflection_points import find_non_singular_inflection_point

try:
    from inflection_points import find_non_singular_inflection_point
except ModuleNotFoundError:
    from .inflection_points import find_non_singular_inflection_point


def multiply(matrix1, matrix2):
    return np.matmul(np.array(matrix1),
                     np.array(matrix2)).tolist() 

def eliminate_denominators(cubic):
    denom = 1
    for coeff in Poly(cubic).coeffs():
        denom = max(denom, abs(coeff.denominator))
    
    return cubic * denom


# In step1 we map inflection point `point` to (0 : 1 : 0)
def weierstrass_form_step1(cubic, point):
    n, x, y, z = symbols('n x y z')
    (xp, yp, zp) = point

    matrix0 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    
    # First, we want to change the variable order, such that if point has zero
    # coordinates, they go like this: (1 : 0 : 0), i.e. if 1 zero, then the
    # point has coordinates (x : y : 0), if 2 zeros, then (x : 0 : 0)
    if zp != 0:
        if xp == 0:
            (zp, xp) = (xp, zp)
            matrix0 = [
                [0, 0, 1],
                [0, 1, 0],
                [1, 0, 0],
            ]

        if yp == 0:
            (zp, yp) = (yp, zp)
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
        (xp, yp) = (yp, xp)
        matrix1 = [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
        ]
    
    
    # Invertable matrix, such that sends our point to (0 : 1 : 0)
    matrix2 = [
        [Fraction(yp, xp), -1, 0],
        [Fraction(1, xp), 0, 0],
        [0, 0, 1],
    ]
    
    matrix = multiply(matrix2,
             multiply(matrix1,
                      matrix0
                      ))
    
    a, b, c = symbols('a b c')

    (x1, y1, z1) = tuple(Matrix(matrix).inv() * Matrix([a, b, c]))

    cubic = cubic.subs({x: x1, y: y1, z: z1}).simplify().expand()
    cubic = cubic.subs({a: x, b: y, c: z})
    
    return (matrix, cubic)



# In step 3 we make a transformation such that z = 0 becomes a tangent to our
# cubic at point (0 : 1 : 0)
def weierstrass_form_step2(cubic):
    n, x, y, z = symbols('n x y z')
    c11, c31, c21, c23, c13, c33 = symbols('c11 c31 c21 c23 c13 c33') 

    # Note, that since after step 1 there is no y^3, then it folows that partial
    # derivative wrt y is always 0
    p = diff(cubic, x).subs({x: 0, y: 1, z: 0})
    q = diff(cubic, z).subs({x: 0, y: 1, z: 0})


    # print(p)
    # print(q)

    if p == 0 and q == 0:
        print("SOMETHING IS WRONG OR UNHANDLED CASE: p == 0 and q == 0")

    
    if p != 0:
        matrix = [
            [-q/p, 0, -q/p - 1],
            [1, 1, 1],
            [1, 0, 1],
        ]


    elif q != 0:
        matrix = [
            [1, 0, 1],
            [1, 1, 1],
            [-p/q, 0, -p/q - 1],
        ]

    # print(matrix)
    # print(diff(cubic, x).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, y).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, z).subs({x: 0, y: 1, z: 0}))

    a, b, c = symbols('a b c')
    (x1, y1, z1) = tuple(Matrix(matrix) * Matrix([a, b, c]))
    
    cubic = simplify(cubic.subs({x: x1, y: y1, z: z1})).expand()
    cubic = cubic.subs({a: x, b: y, c: z})

    (matrix1, cubic) = simplify_cubic(cubic)
    matrix = Matrix(matrix).inv().tolist()
    

    # print(diff(cubic, x))
    # print(diff(cubic, x).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, y).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, z).subs({x: 0, y: 1, z: 0}))

    return (multiply(matrix1, matrix), cubic)


# This function just tries to reduce coefficients as much as possible
# Matrix is the actual matrix in projective transformation (not inverted)
#
# TODO: maybe I will need to involve z in simplification
def simplify_cubic(cubic):
    n, x, y, z = symbols('n x y z')
    z_, y_ = symbols('z_ y_')
    
    cubic = eliminate_denominators(cubic)

    a = abs(cubic.coeff(y**2 * z))
    primes = factorint(a)
    coeff = 1

    for i in primes.keys():
        if primes[i] % 2 == 0:
            coeff *= i**(int(primes[i]/2))
    
    cubic = cubic.subs(y, y_/coeff).simplify().expand()
    cubic = cubic.subs(z, z_ * coeff**2/a).simplify().expand()
    
    matrix = [
        [1, 0, 0],
        [0, coeff, 0],
        [0, 0, a/coeff**2],
    ]

    cubic = cubic.subs(y_, y).simplify().expand()
    cubic = cubic.subs(z_, z).simplify().expand()

    return (matrix, cubic)


def to_integer_matrix(matrix):
    mat_list = np.array(matrix).flatten().tolist()
    
    lcm = np.lcm.reduce([i.denominator for i in mat_list])
    mat_list = (Matrix(mat_list) * lcm).tolist()
    gcd = np.gcd.reduce([i for i in mat_list])

    matrix = Matrix(matrix) * lcm/gcd
    
    return matrix.tolist()


def weierstrass_form_step3(cubic):
    n, x, y, z = symbols('n x y z')
    x_, y_, z_ = symbols('x_ y_ z_')

    coeff = cubic.coeff(y**2 * z)
    cubic = cubic / coeff

    a = 1
    b = cubic.coeff(y * x * z)
    c = cubic.coeff(y * z**2)
    
    # h, k, t = symbols('h k t')
    # (p, q)  = tuple(solve((t + h)**2 + k  - (t**2  + b * t * x  + c * t * z), [h, k])[0])
    p = (b * x + c * z)/2
    # print(p, q)

    # print(solve((t + h)**2 + k  - (t**2  + b * t * x  + c * t * z), [h, k]))
    # print(cubic)

    cubic = cubic.subs(y, y_ - p).expand().simplify().expand()
    cubic = cubic.subs(y_, y)
    cubic = (-1) * cubic
    
    # Complete the square by `y`
    matrix0 = [
        [1, 0, 0],
        [p.coeff(x), 1, p.coeff(z)],
        [0, 0, 1],
    ]

    # Reduce coefficients, because otherwise they are enormous
    (matrix1, cubic) = simplify_cubic(cubic)

    
    a = cubic.coeff(x**3)
    b = cubic.coeff(x**2 * z)
    # c = cubic.coeff(x * z**2)
    # d = cubic.coeff(z**3)

    if (a == 0):
        print("UNHANDLED CASE a == 0, step 3")
    
    # print(Fraction(b, 3 * a))

    
    # Eliminating the x^2 z monomial
    cubic = cubic.subs(x, x_ - Fraction(b, 3 * a) * z).expand().simplify().expand()
    cubic = cubic.subs(x_, x)
    matrix2 = [
        [1, 0, Fraction(b, 3 * a)],
        [0, 1, 0],
        [0, 0, 1],
    ]
    

    # Normalizing the coefficient of x^3
    a = cubic.coeff(x**3)
    cubic = (cubic.subs(z, z_ * a) / a).expand().simplify().expand()
    cubic = cubic.subs(z_, z)
    matrix3 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, Fraction(1, a)],
    ]
   
    
    lcm = np.lcm.reduce([coeff.denominator for coeff in Poly(cubic).coeffs()])

    # Getting rid of denominators in the fractions
    cubic = cubic.subs(z, z_ * lcm**3).expand().simplify().expand()
    cubic = cubic.subs(x, x_ * lcm).expand().simplify().expand()
    matrix4 = [
        [Fraction(1, lcm), 0, 0],
        [0, 1, 0],
        [0, 0, Fraction(1, lcm**3)],
    ]
    
    matrix = multiply(matrix4,
             multiply(matrix3,
             multiply(matrix2,
             multiply(matrix1,
                      matrix0
                      ))))

    cubic = (cubic / lcm**3).simplify().expand()
    cubic = cubic.subs({x_: x, z_: z})

    return (matrix, cubic)


    
def weierstrass_form(cubic):
    n, x, y, z = symbols('n x y z')

    (res, point) = find_non_singular_inflection_point(cubic)
    if not res:
        return (False, []) 
        
    (trans1, cubic) = weierstrass_form_step1(cubic, point)
    (trans2, cubic) = weierstrass_form_step2(cubic)
    (trans3, cubic) = weierstrass_form_step3(cubic)
    # print(cubic)
    # return 

    trans = multiply(trans3,
            multiply(trans2,
                     trans1
                     ))

    a = cubic.coeff(x * z**2)
    b = cubic.coeff(z**3)
    
    trans = to_integer_matrix(trans)

    return (True, ((a, b), trans))


def print_form(form):
    ((a, b), trans) = form

    a_str = ""
    b_str = ""

    if b < 0:
        b_str = " - {} z^3".format(-b)
    
    if b > 0:
        b_str = " + {} z^3".format(b)
    
    if a < 0:
        a_str = " - {} x z^2".format(-a)

    if a > 0:
        a_str = " + {} x z^2".format(a)


    print("""Cubic equation in Weierstrass form:\n
        y^2 z = x^3{}{}
    """.format(a_str, b_str))

    print("""Corresponding transformation matrix:\n """)

    pprint(Matrix(trans))

    print("\n")


def main():
    print("""Enter your cubic's equation in homogenious coordinates x, y, z:
    For example: x^3 + y^3 + z^3 + 3 x y z
    """)

    n, x, y, z = symbols('n x y z')

    cubic = input()
    cubic = mathematica(cubic)

    (res, form) = weierstrass_form(cubic)
    
    if res:
        print_form(form)
    else:
        print("Couldn't apply algorithm due to previous error")



if __name__ == '__main__':
    main()

