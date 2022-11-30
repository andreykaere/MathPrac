#!/usr/bin/python3

# Finding parent diretory
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ..weierstrass_form.weierstrass_main import weierstrass_form


from sympy import *
from sympy.parsing.mathematica import mathematica
from sympy.ntheory import factorint
from fractions import Fraction
import numpy as np

# import importlib
# importlib.import_module('..inflection_points')
# from inflection_points import find_non_singular_inflection_point



# Takes cubic as a string written in Mathematica style
def apply_test(cubic):
    cubic = mathematica(cubic)
   
    (res, trans) = find_non_singular_inflection_point(cubic)

    if res:
        pprint(Marix(trans))
    else:
        print("Couldn't apply algorithm to the given cubic")


def main_test():
    cubic = "x^3 + y^3 + z^3 + (1 - n) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 n) x y z"

def test1():
    cubic = "-x^3 - x^2*z + y^2*z + 2*y*z^2 + z^3"
    return cubic

def test2():
    cubic = "y^2 z - x^3 - x^2 z"
    return cubic
    
def test3():
    cubic = "5 y^3 + z^2 x + y^2 x - 34 y^2 z"
    return cubic

def test4():
    cubic = "(x - z) (x z - y^2)"
    return cubic

def test5():
    cubic = "(x - y) (y^2 - x^2 + z x) - x^2 y"
    return cubic

def test6():
    cubic = "x^3 + y^3 + z^3 + 3 x y z"
    return cubic

def test7():
    cubic = "x^3 - y^2 z"
    return cubic

def test8():
    cubic = "(z + y) (3 x - y) (x + 3 z)"
    return cubic

def test9():
    cubic = 


def main():
    n, x, y, z = symbols('n x y z')
    
    apply_test(test1())





if __name__ == "__main__":
    main()
