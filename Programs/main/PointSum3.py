from fractions import Fraction
import matplotlib.pyplot as plt
import math
from math import *
import numpy as np
import itertools
import matplotlib.animation as animation
from time import *



# only for z = 1
# int a, b ;      Fraction xp, yp, xq, yq  better for accuracy 
# but mb double
def PointSum(a, b, xp, yp, xq, yq):
    # https://habr.com/ru/post/335906/?ysclid=l9bh1vdtjn923540758
    if (yp**2 != xp**3 + a * xp + b):
        print("!!!!Point not on curve!!!!")
        return -1
    if (yq**2 != xq**3 + a * xq + b):
        print("!!!!Point not on curve!!!!")
        return -1
    #dont work when tangent is verlical
    if (xp == xq) and (yp == yq):
        if (yp == 0):
            print("!!Tangent is vertical!!")
        m = (3 * xp * xp + a) / (2 * yp)
        xr = m * m - xp - xq
        yr = yp + m * (xr - xp)           

    else:
        m = (yp - yq)/ (xp - xq)
        xr = m * m - xp - xq
        yr = yp + m * (xr - xp)   
    
    return [xr, -yr]

#for P = (x, y) return n*P = P + P +... + P
def ScMult(a, b, x, y, n):
    x0 = x
    y0 = y
    #print(1,x, y)
    for i in range(n - 1):
        A = PointSum(a, b, x, y, x0, y0)
        x = A[0]
        y = A[1]
        #print(i + 2, x, y)
    return [x, y]

#for point P = (x, y) return list of n point: P, P+P, .. n*P
def GenerateNpoints(a, b, x, y, n):
    xP = [x]
    yP = [y]
    x0 = x
    y0 = y
    for i in range(2, n + 1):
        x0 = ScMult(a, b, x, y, i)[0]
        y0 = ScMult(a, b, x, y, i)[1]
        xP.append(x0)
        yP.append(y0)
    return [xP, yP]


#https://stackoverflow.com/questions/2489435/check-if-a-number-is-a-perfect-square
def is_square(apositiveint):
    if (apositiveint == 1):
        return True
    x = apositiveint // 2
    seen = set([x])
    while x * x != apositiveint:
        x = (x + (apositiveint // x)) // 2
        if x in seen: return False
        seen.add(x)
    return True

# Finds (if possible) in range R for numerator and denominator
def FindRational(a, b, R):
    for xN in range(-R, R):
        for xD in range(1, R):
            if gcd(abs(xN), xD) == 1:
                x = Fraction(xN, xD)
                y2 = (x**3 + a * x + b)
            
                y2n =  abs(y2.numerator)
                y2d = abs(y2.denominator)
                if (y2n == 0):
                    return [[x, Fraction(0)], 1]
                    
                if (y2n != 0) and is_square(y2n) and is_square(y2d) and (y2 > 0):
                    return [[x, Fraction(int(y2n**(1/2)), int(y2d**(1/2)))], 1]
                    
    print("END")
    return [[], 0]


#finds coordinate of point in the original coordinate system
def Reverse(x, y, z, N):
    N = Fraction(N)
    #to do: chech if point on curve
    a = [[1/(72*(3 + N)*(5 + 2*N)), 
     1/(216*(3 + N)*(5 + 2*N)), 
     -((69 + 4*N*(9 + N))/(24*(3 + N)*(5 + 2*N)))], 
   [1/(72*(3 + N)*(5 + 2*N)), 
     -(1/(216*(3 + N)*(5 + 2*N))), 
     -((69 + 4*N*(9 + N))/(24*(3 + N)*(5 + 2*N)))], 
   [1/(36*(7 + 2*N + 1/(2 + N))), 0, 
     Fraction(1,12)*(1 - 2*N + 3/(3 + N) + 4/(5 + 2*N))]]
    
    x0 = a[0][0] * x + a[0][1] * y + a[0][2] * z
    y0 = a[1][0] * x + a[1][1] * y + a[1][2] * z
    z0 = a[2][0] * x + a[2][1] * y + a[2][2] * z
    return [x0, y0, z0]

# finds first point from P_1, .. P_n,
# that has natural coordinates in the original coordinate system
# and return these coordinates
def RevFind(xP, yP, N, start = 0):
    n = len(xP)
    zP = []
    for i in range(n):
        zP.append(Fraction(1))

    for i in range(start, n):
        x = xP[i]
        y = yP[i]
        z = zP[i]
        V = Reverse(x, y, z, N)
        x0 = V[0]
        y0 = V[1]
        z0 = V[2]
    
        if (x0 >= 0) and (y0 >= 0) and (z0 >= 0):

            return [[x0, y0, z0], i]
        if (x0 <= 0) and (y0 <= 0) and  (z0 <= 0):
            return [[x0, y0, z0], i]
    print("No more points")
    return 0

# for point P = (x, y) finds 2*P, .. n*P
# and finds all points that have natural coordinate
# in the original system 
# return their coordinates
# it is the solution for our equation
def RevNFind(a, b, x, y, n, N):
    
    V = GenerateNpoints(a, b, x, y, n)
    xP = V[0]
    yP = V[1]
    PointNum = 0
    start = 0
    while True:
        Ans = RevFind(xP, yP, N, start)
        if type(Ans) == int:
            break
        PointNum += 1
        a = Ans[0][0]
        b = Ans[0][1]
        c = Ans[0][2]
        start = Ans[1] + 1

        ad = a.denominator
        bd = b.denominator
        cd = c.denominator
        sgn = 1
        if (a < 0):
            sgn = -1
        a = a.numerator * bd * cd * sgn
        b = b.numerator * ad * cd * sgn
        c = c.numerator * ad * bd * sgn
        d = gcd(gcd(a,b), c)
        a = a // d
        b = b // d
        c = c // d
        if PointNum == 1:
            minLen = min(len(str(a)), len(str(b)), len(str(c)))
        else:
            minLen = min(len(str(a)), len(str(b)), len(str(c)), minLen)

  
        print(PointNum, "------------------------------------")
        print("a =", a)
        print()
        print("b =", b)
        print()
        print("c =", c)

        a = Fraction(a, 1)
        b = Fraction(b, 1)
        c = Fraction(c, 1)

        print()
        print("a/(b+c) + b/(a+c) + c/(a+b) =", a/ (b + c) + b/(a + c) + c/ (a + b))
        print("-------------------------------")
   
    print()
    print("Points found:", PointNum)
    print("Min length", minLen)


##-----------------------
#N = 4
#a = -432 * (N**4) - 2592 * (N**3) - 3240 * (N**2) + 4536 * N + 7533
#b = 3456 * (N**6) + 31104 * (N**5) + 85536 * (N**4) + 15552 * (N**3) - 250776 * (N**2) - 239112 * N + 68526
#x = Fraction(246, 1)
#y = Fraction(2106, 1)


#n = 9



#DrawSumN(a, b, x, y, n)

##RevNFind(a, b, x, y, n)
##-----------------------

