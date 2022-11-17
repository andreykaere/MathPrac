#!/usr/bin/python3

from sympy import *
# from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint
# from sympy.ntheory import primefactors

from fractions import Fraction
import numpy as np


def find_common_projective_solutions(f, g):
    n, x, y, z = symbols('n x y z')
    t = symbols('t')

    solutions = set()


    
    # If x != 0
    res1 = resultant(f, g, z).simplify().expand()
    k1 = Poly(res1).homogeneous_order()
    res_t1 = collect(expand(res1/x**k1).subs(y/x, t), t).simplify().expand()

    print(solve(res_t1, t))
    print(res_t1)


    # If y != 0
    res_t2 = collect(expand(res1/y**k1).subs(x/y, t), t).simplify().expand()

    print(solve(res_t2, t))
    print(res_t2)



    

# TODO: fix it 
def is_cubic_singular(cubic):
    return False 

    n, x, y, z = symbols('n x y z')
    t = symbols('t')
    
    fx = diff(cubic, x)
    fy = diff(cubic, y)
    fz = diff(cubic, z)

    points1 = set(find_common_projective_solutions(fx, fy))
    points2 = set(find_common_projective_solutions(fx, fz))
    singular_points = list(points1.intersection(point2))

    if singular_points != []:
        print("Cubic is singular at points", singular_points)
        return True

    return False


    # First, check if there is a singular point (x0 : y0 : 0) 
    # fx0 = fx.subs(z, 0).simplify().expand()
    # fy0 = fy.subs(z, 0).simplify().expand()
    # fz0 = fz.subs(z, 0).simplify().expand()
    # print(fx0)
    # print(fy0)
    # print(fz0)



    # g0 = resultant(fx0, fy0, x).simplify().expand()
    # h0 = resultant(fx0, fz0, x).simplify().expand()
    # k0 = Poly(h0).homogeneous_order()
    # m0 = Poly(g0).homogeneous_order()

    # print(g0)
    # print(h0)

    # return
    
    # Now x != 0 or y != 0
    # First, if y != 0 
    # h0_xy_t = collect(expand(h0/y**k).subs(x/y, t), t).simplify().expand()
    # g0_xy_t = collect(expand(g0/y**m).subs(x/y, t), t).simplify().expand()
    # res0_xy = resultant(h0_xy_t, g0_xy_t, t).simplify().expand()
    
    # # Now if x != 0
    # h0_yx_t = collect(expand(h0/y**k).subs(y/x, t), t).simplify().expand()
    # g0_yx_t = collect(expand(g0/y**m).subs(y/x, t), t).simplify().expand()
    # res0_yx = resultant(h0_yx_t, g0_yx_t, t).simplify().expand()

    # print(res0_xy)
    # print(res0_yx)

    # if cubic.subs({x: 0, y: 0, z: 1}) == 0:
    #     # print("point (0 : 0 : 1) is singular")
    #     return True

#     g = resultant(fx, fy, z).simplify().expand()
#     h = resultant(fx, fz, z).simplify().expand()
#     k = Poly(h).homogeneous_order()
#     m = Poly(g).homogeneous_order()
    
#     # Now x != 0 or y != 0
#     # First, if y != 0 
#     h_xy_t = collect(expand(h/y**k).subs(x/y, t), t).simplify().expand()
#     g_xy_t = collect(expand(g/y**m).subs(x/y, t), t).simplify().expand()
#     res_xy = resultant(h_xy_t, g_xy_t, t).simplify().expand()
    
#     # Now if x != 0
#     h_yx_t = collect(expand(h/x**k).subs(y/x, t), t).simplify().expand()
#     g_yx_t = collect(expand(g/x**m).subs(y/x, t), t).simplify().expand()
#     res_yx = resultant(h_yx_t, g_yx_t, t).simplify().expand()
    
    # print(g)
    # print(h)
    # print(k, m)
    # print(h_xy_t)
    # print(g_xy_t)
    # print(h_yx_t)
    # print(g_yx_t)
    # print(res_xy, res_yx)

    # If one of these is zero, than there is a common solutions for fx, fy,
    # fz: either (0, 1) or (1, 0) therefore there a point (0 : 1 : z0) or 
    # (1 : 0 : z0) such that all fx, fy, fz is zero, i.e. singular point
    # if res_xy == 0 or res_yx == 0:
    #     return True
    
    # return False

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


def from_solution_to_rational_points(cubic, sol):
    n, x, y, z, t = symbols('n x y z t')
    # x_, y_, z_ = symbols('n x y z')
    xp = sol[0]
    yp = sol[1]

    cubic = cubic.subs({x: xp, y: yp}).simplify().expand()
    rational_z = get_rational_roots(Poly(cubic).subs(z, t))

    if rational_z == []:
        return (False, [])    
    
    return (True, [(xp, yp, zp) for zp in rational_z if (xp, yp, zp) != (0, 0, 0)])
     
    
def fix_zero_leading_coefficients(cubic, hessian):
    n, x, y, z = symbols('n x y z')
    i = 1 

    # Searching for point that does not lie on both cubics
    while cubic.subs({x: 0, y: i, z: 1}) == 0 and hessian.subs({x: 0, y: i, z: 1}) == 0:
        i += 1
   
    point = (0, i, 1)
    
    matrix = [
        [1, 0, 0],
        [0, 1, i],
        [1, 0, 1],
    ]

    
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
    coeff_t = poly_t.all_coeffs()
    solutions = set()

    if poly_t.subs(t, 0) == 0:
        solutions.add(0)

    # Getting rid of multiple zero roots so that we can start our enumeration 
    # algorithm 
    while poly_t.subs(t, 0) == 0:
        poly_t = expand(poly_t / t)

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
    

def points_of_inflection(cubic):
    n, x, y, z = symbols('n x y z')
    trans = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    hessian = get_hessian(cubic)

    a0 = cubic.coeff(z**3)
    b0 = hessian.coeff(z**3)

    if a0 * b0 == 0:
        (cubic, hessian, trans) = fix_zero_leading_coefficients(cubic, hessian)

    t = symbols('t')
    res = resultant(cubic, hessian, z)
    
    # Creating set and not array, because we don't care if roots are multiple 
    # or not
    solutions = {(0, 0)}

    if res.subs({x: 1, y: 0}) == 0:
        solutions.add((1, 0))


    res_t = Poly(collect(expand(res/y**9).subs(x/y, t), t), t)
    # print(res_t)
    
    # Reducing our polynom, so that our enumeration algorithm will not take 
    # forever to finish
    coeff_t = res_t.all_coeffs()
    gcd = np.gcd.reduce(coeff_t)
    res_t = simplify(res_t / gcd).expand()


    if res_t.subs(t, 0) == 0:
        solutions.add((0, 1))
        res_t = expand(res_t / t)
    
    
    rational_sols = get_rational_roots(Poly(res_t))

    for sol in rational_sols:
        solutions.add((sol, 1))

    
    points = []
    
    for solution in list(solutions):
        (result, points_loc) = from_solution_to_rational_points(cubic, solution)

        # print(points_loc)
        
        if result:
            points += points_loc
    
    if points == []:
        print("Fatal: No rational inflection points exist on the cubic, can't proceed")

    # print(points)

    return (points, trans)


def main():
    n, x, y, z = symbols('n x y z')

    # cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"
    # cubic = "y^2 z - x^3 - x^2 z"
    # cubic = "-x^3 - x^2*z + y^2*z + 2*y*z^2 + z^3"
    # cubic = "-x^3 - 4*x^2*z + y^2*z - 5*x*z^2 - 2*z^3"
    # cubic = "-x^3 - 4*x^2*z + y^2*z - 5*x*z^2 + 2*y*z^2 - z^3"
    # cubic = "-x^3 - 3*x^2*z + y^2*z - 3*x*z^2"
    cubic = "x^3 + y^2 z + z^3"
    cubic = mathematica(cubic)
    # cubic = cubic.subs(n, 4)

    print(is_cubic_singular(cubic))
    
    # print(cubic)
    # print(get_hessian(cubic))

    # print(points_of_inflection(cubic))



if __name__ == "__main__":
    main()
