#!/usr/bin/python3

import time
from math import sqrt


def is_square(x):
    for i in range(int(sqrt(x)) + 1):
        if i**2 == x:
            return True
    return False


def is_square_mod(x, a):
    for i in range(a):
        if (i**2 - x) % a == 0:
            return True
    return False


def CL(a, b, c):
    return is_square_mod(b * c, a) and \
           is_square_mod(a * c, b) and \
           is_square_mod(-a * b, c)


def find_solution(a, b, c):
    x_range = int((b * c) ** (1/2)) + 1
    y_range = int((a * c) ** (1/2)) + 1
    z_range = int((a * b) ** (1/2)) + 1

    # mb not enough precision
    for x0 in range(x_range):
        for y0 in range(y_range):
            z0 = (a * x0**2 + b * y0**2)/c

            if is_square(z0) and (x0 != 0 or y0 != 0 or z0 != 0):
                return (x0, y0, z0)

    
    

def main():
    print("input a, b, c")

    a = int(input())
    b = int(input())
    c = int(input())

    start_time = time.time()
    if (not CL(a, b, c)):
        print("No solution")
    
    else:
        print("------ %s seconds ------" % (time.time() - start_time))
        print("There is a solution, wait....................")
        
        print(find_solution(a, b, c))
    
    print("-------- %s seconds -------" % (time.time() - start_time))
    print()


main()
