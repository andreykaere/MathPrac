#!/usr/bin/python3

from sympy import *
# from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint
# from sympy.ntheory import primefactors

from fractions import Fraction
import numpy as np

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

    matrix2 = [
        [Fraction(yp, xp), -1, 1],
        [Fraction(1, xp), 0, 0],
        [0, 0, 1],
    ]
    
    # matrix = np.dot(matrix2,
    #          np.dot(matrix1, 
    #                 matrix0
    #          ).tolist()
    #          ).tolist()
    
    matrix = multiply(matrix2,
             multiply(matrix1,
                      matrix0
                      ))
    
    a, b, c = symbols('a b c')

    (x1, y1, z1) = tuple(Matrix(matrix).inv() * Matrix([a, b, c]))

    cubic = simplify(cubic.subs({x: x1, y: y1, z: z1}))
    cubic = cubic.subs({a: x, b: y, c: z})
    
    # print(cubic)

    return (matrix, cubic)

    

def weirstrass_form_step2(cubic):
    n, x, y, z = symbols('n x y z')
    c11, c31, c21, c23, c13, c33 = symbols('c11 c31 c21 c23 c13 c33') 

    # matrix = Matrix([
    #     [c11, 0, c13],
    #     [c21, 1, c23],
    #     [c31, 0, c33],
    # ])
    

    # Note, that since after step 1 there is no y^3, then it folows that partial
    # derivative wrt y is always 0
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
        print("SOMETHING IS WRONG OR UNHANDLED CASE: p == 0 and q == 0")

    
    if p != 0:
        matrix = [
            [-q/p, 0, -q/p - 1],
            [1, 1, 1],
            [1, 0, 1],
        ]
        # matrix = [
        #     [-q, 0, -q - 1],
        #     [1, 1, 1],
        #     [p, 0, p],
        # ]


    elif q != 0:
        matrix = [
            [1, 0, 1],
            [1, 1, 1],
            [-p/q, 0, -p/q - 1],
        ]
        # matrix = [
        #     [-q, 0, -q],
        #     [1, 1, 1],
        #     [p, 0, p + 1],
        # ]



    # print(matrix)
    # print(diff(cubic, x).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, y).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, z).subs({x: 0, y: 1, z: 0}))

    a, b, c = symbols('a b c')
    (x1, y1, z1) = tuple(Matrix(matrix) * Matrix([a, b, c]))
    
    cubic = simplify(cubic.subs({x: x1, y: y1, z: z1}))
    cubic = cubic.subs({a: x, b: y, c: z})

    
    # print(matrix)
    # print(Matrix(matrix).inv().tolist())
    # print(cubic)


    # print(diff(cubic, x))
    # print(diff(cubic, x).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, y).subs({x: 0, y: 1, z: 0}))
    # print(diff(cubic, z).subs({x: 0, y: 1, z: 0}))

    return (Matrix(matrix).inv().tolist(), cubic)




def eliminate_denominators(cubic):
    denom = 1

    for coeff in Poly(cubic).coeffs():
        denom = max(denom, abs(coeff.denominator))
    
    return cubic * denom




# This function just tries to reduce coefficients as much as possible
# //We also assume here, that `eliminate_denominators` function was called if
# //needed
#
# TODO: maybe I will need to involve z in simplification
def simplify_weirstrass_form(cubic):
    n, x, y, z = symbols('n x y z')
    z_, y_ = symbols('z_ y_')
    
    cubic = eliminate_denominators(cubic)

    a = abs(cubic.coeff(y**2 * z))

    # print(a)

    # a = decompose(a) ..
    primes = factorint(a)

    coeff = 1

    for i in primes.keys():
        if primes[i] % 2 == 0:
            coeff *= i**(int(primes[i]/2))
    
    cubic = cubic.subs(y, y_/coeff).simplify()
    cubic = cubic.subs(z, z_ * coeff**2/a).simplify()
    
    matrix = [
        [1, 0, 0],
        [0, coeff, 0],
        [0, 0, a/coeff**2],
    ]

    cubic = cubic.subs(y_, y).simplify()
    cubic = cubic.subs(z_, z).simplify()

    return (matrix, cubic)


def multiply(matrix1, matrix2):
    return np.dot(matrix1, matrix2).tolist() 

def simplify_transformation(matrix):
    # print(matrix) 
    mat_list = np.array(matrix).flatten().tolist()

    # print(mat_list)

    gcd = np.gcd.reduce([i.denominator for i in mat_list])
    # print([i.denominator for i in mat_list])

    print(Matrix(mat_list) * gcd)

    return matrix

def weirstrass_form_step3(cubic):
    n, x, y, z = symbols('n x y z')
    x_, y_, z_ = symbols('x_ y_ z_')

    coeff = cubic.coeff(y**2 * z)
    cubic = cubic / coeff

    a = 1
    b = cubic.coeff(y * x * z)
    c = cubic.coeff(y * z**2)
    
    h, k, t = symbols('h k t')
    (p, q)  = tuple(solve((t + h)**2 + k  - (t**2  + b * t * x  + c * t * z), [h, k])[0])
    # print(p, q)

    # print(solve((t + h)**2 + k  - (t**2  + b * t * x  + c * t * z), [h, k]))
    # print(cubic)

    cubic = cubic.subs(y, y_ - p).expand().simplify()
    cubic = cubic.subs(y_, y)
    
    # Complete square by `y`
    matrix0 = [
        [1, 0, 0],
        [p.coeff(x), 1, p.coeff(z)],
        [0, 0, 1],
    ]


    # print(cubic)
    # print(matrix0, cubic)
    (matrix1, cubic) = simplify_weirstrass_form(cubic)

    print(matrix1, cubic)



    
    a = cubic.coeff(x**3)
    b = cubic.coeff(x**2 * z)
    # c = cubic.coeff(x * z**2)
    # d = cubic.coeff(z**3)

    if (a == 0):
        print("UNHANDLED CASE a == 0, step 3")
    
    # print(Fraction(b, 3 * a))


    cubic = cubic.subs(x, x_ - Fraction(b, 3 * a) * z).expand().simplify()
    cubic = cubic.subs(x_, x)
    matrix2 = [
        [1, 0, b/(3 * a)],
        [0, 1, 0],
        [0, 0, 1],
    ]
    

   
    a = cubic.coeff(x**3)
    cubic = (cubic.subs(z, z_ * a) / a).expand().simplify()
    cubic = cubic.subs(z_, z)
    matrix3 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1/a],
    ]
    # return (matrix3, cubic)
    # print(cubic)

    
    denom = 1

    for coeff in Poly(cubic).coeffs():
        denom = max(denom, abs(coeff.denominator))

    cubic = cubic.subs(z, - z_ * denom**3).expand().simplify()
    cubic = cubic.subs(x, x_ * denom).expand().simplify()
    matrix4 = [
        [Fraction(1, denom), 0, 0],
        [0, 1, 0],
        [0, 0, -Fraction(1, denom**3)],
    ]
    
    matrix = multiply(matrix4,
             multiply(matrix3,
             multiply(matrix2,
             multiply(matrix1,
                      matrix0
                      ))))

    # matrix = simplify_transformation(matrix)

    cubic = (cubic / denom**3).simplify()

    cubic = cubic.subs({x_: x, z_: z})
    # print(cubic)

    return (matrix, cubic)



    
def weirstrass_form(cubic):
    n, x, y, z = symbols('n x y z')

    # TODO: uncomment
    # (points, trans0) = points_of_inflection(cubic)
    # point = points[0]
    
    point = (1, -1, 0)
    trans0 = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    
    (trans1, cubic) = weirstrass_form_step1(cubic, point)
    (trans2, cubic) = weirstrass_form_step2(cubic)
    (trans3, cubic) = weirstrass_form_step3(cubic)

    # print(trans1)
    # print(trans2)
    # print(trans3)
   
    # trans = np.dot(trans3,
    #         np.dot(trans2, 
    #         np.dot(trans1,
    #                trans0
    #         ).tolist()
    #         ).tolist()
    #         ).tolist() 

    trans = multiply(trans3,
            multiply(trans2,
            multiply(trans1,
                     trans0
                     )))

    a = cubic.coeff(x * z**2)
    b = cubic.coeff(z**3)

    # print(cubic)
    # print(a,b)

    return ((a, b), trans)

    





def main():
    # print("Enter your cubic's equation in homogenious coordinates x, y, z:")

    # cubic = input()
    n, x, y, z = symbols('n x y z')

    cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"

    cubic = mathematica(cubic)

    cubic = cubic.subs(n, 4)

    print(weirstrass_form(cubic))

    # print(weirstrass_form_step1(cubic, (-1, 1, 0)))
    # (trans1, cubic) = weirstrass_form_step1(cubic, (-1, 1, 0))
    # print(trans1, cubic)
    # (trans2, cubic) = weirstrass_form_step2(cubic)
    # (trans3, cubic) = weirstrass_form_step3(cubic)
    # print(trans3, cubic)
    # print(cubic)
    # print(weirstrass_form_step3(cubic))


if __name__ == '__main__':
    main()

# print(resultant(get_gessian(cubic), cubic, z))
