#!/usr/bin/python3

from fractions import Fraction
import matplotlib.pyplot as plt
import math
from math import *
import numpy as np
import itertools
import matplotlib.animation as animation


# only for z = 1, for now
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
        #print("__________________________")
        #print(xp, yp)
        #print(xq, yq)
        #print("__________________________")
        #print()
        m = (yp - yq)/ (xp - xq)
        xr = m * m - xp - xq
        yr = yp + m * (xr - xp)   
    
    return [xr, -yr]

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

#return list of n point: P, P+P, .. n*P
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
    
def GenerateAlotPoints(a, b, x, y, n = 20):
    xP = [x]
    yP = [y]
    x0 = x
    y0 = y    
    for i in range(n - 1):
        try:
            A = PointSum(a, b, x, y, x0, y0)
        except:
            A = PointSum(a, b, x, y, xP[1], yP[1])
            
        x = A[0]
        y = A[1]
        xP.append(x)
        yP.append(y)
    return [xP, yP]  

def Draw(a = -5, b = 5, xPoints = [], yPoints = []):
    x = np.arange(-10,10,0.01)
    yu = ((x**3 + a * x + b))**(1/2)
    yd = - ((x**3 + a * x + b))**(1/2)
    plt.title('My first plot')
    plt.plot(x,yu,'b', x, yd, 'g', xPoints, yPoints, 'ro')
    plt.show()    

# to draw plt.show()
def DrawCurve(a, b):
    x = np.arange(-1000,1000,0.1)
    yCurveU = []
    yCurveD = []
     
    for i in x.tolist():
        if i**3 + a*i + b >= 1/100:    
            yCurveU += [((i**3 + a * i + b))**(1/2)]
            yCurveD += [-((i**3 + a * i + b))**(1/2)]
        else:
            i0 = int(x.tolist().index(i))
            #x = np.delete(x, i0, 0)
            x[i0] = np.nan
            yCurveU += [np.nan]
            yCurveD += [np.nan]            
            #yCurveU += [0]
            #yCurveD += [0]
   
    yCurveU = np.array(yCurveU)
    yCurveD = np.array(yCurveD)    
    
    plt.plot(x,yCurveU,'b', x, yCurveD, 'b')
    
def DrawSum(a = -5, b = 5, xPoints = [], yPoints = []):

    xL = np.arange(-1000,1000,0.1)
 
    xp = xPoints[0]
    yp = yPoints[0]
    xq = xPoints[1]
    yq = yPoints[1] 
    
    A = PointSum(a, b, xp, yp, xq, yq)

    xr = A[0]
    yr = A[1]
    xr = [xr]
    yr = [yr]
    
    if (xp == xq) and (yp == yq):
        m = (3 * xp * xp + a) / (2 * yp)
    
    else:
        m = (yp - yq)/ (xp - xq)
    yLine = m * (xL - xq) + yq
    xLine1 = [xr]*60000
    yLine1 = np.arange(-30000, 30000, 1)    
    plt.title('Sum of points')

    plt.plot(xL, yLine, linestyle = '-', linewidth = 1, color = 'y')
    plt.plot(xLine1, yLine1, linestyle = '--', linewidth = 1, color = 'darkmagenta')
    plt.plot(xPoints, yPoints, 'ro', xr, yr, 'go')
    #plt.show()

# (not used) for x find y(x) > 0 
def yC(x, a, b):
    y2 = (x**3 + a * x + b)
    y2n =  abs(y2.numerator)
    y2d = abs(y2.denominator)
    if (abs(y2n) > 10**8) or (abs(y2d) > 10**8):
        print("!!!Check square mb unsave!!!")
    
    if (int(y2n**(1/2)) == y2n**(1/2)) and (int(y2d**(1/2)) == y2d**(1/2)):
            return Fraction(y2**(1/2))

    else:
        print("!! Not square!!")
        return y2**(1/2)

def DrawSumN(a, b, x, y, n):
    xP = [x]
    yP = [y]
    x0 = x
    y0 = y
    fig, ax = plt.subplots(figsize=(8, 8))
    
    DrawCurve(a, b)
    for i in range(1, n + 1):
        DrawSum(a, b, [x, x0], [y, y0])
        if (i != 1):
            plt.text(x0, y0 + 1, str(i), fontsize=9)
        x0 = ScMult(a, b, x, y, i)[0]
        y0 = ScMult(a, b, x, y, i)[1]
        xP.append(x0)
        yP.append(y0)
    plt.plot([x], [y], 'yo')
    #plt.text(x0, y0+1, "result", fontsize=9)
    plt.show() 
#https://stackoverflow.com/questions/2489435/check-if-a-number-is-a-perfect-square
def is_square(apositiveint):
    if (apositiveint == 1):
        return True
    x = apositiveint // 2
    seen = set([x])
    while x * x != apositiveint:
        #print(x)
        try:
            x = (x + (apositiveint // x)) // 2
        except:
            print(x)
        if x in seen: return False
        seen.add(x)
    return True

def FindRational(a, b, R):
    for xN in range(0, R):
        for xD in range(1, R):
            x = Fraction(xN, xD)
            y2 = (x**3 + a * x + b)
            
            y2n =  abs(y2.numerator)
            y2d = abs(y2.denominator)
            if (y2n == 0):
                print([x, Fraction(0)])
            if (y2n != 0) and is_square(y2n) and is_square(y2d):
                print( [x, Fraction(int(y2n**(1/2)), int(y2d**(1/2)))])
    print("END")

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

def RevFind(xP, yP, N):
    n = len(xP)
    zP = []
    for i in range(n):
        zP.append(Fraction(1))

    for i in range(n):
        x = xP[i]
        y = yP[i]
        z = zP[i]
        V = Reverse(x, y, z, N)
        x0 = V[0]
        y0 = V[1]
        z0 = V[2]
    
        if (x0 >= 0) and (y0 >= 0) and (z0 >= 0):
            return [x0, y0, z0]
        if (x0 <= 0) and (y0 <= 0) and  (z0 <= 0):
            return [x0, y0, z0]
    print("No points")
    return 0

<<<<<<< HEAD
def RevNFind(a, b, x, y, n):
    
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
=======

>>>>>>> parent of 111ee21 (added some things to SumPoint)
#task data
N = 4
a = -432 * (N**4) - 2592 * (N**3) - 3240 * (N**2) + 4536 * N + 7533
b = 3456 * (N**6) + 31104 * (N**5) + 85536 * (N**4) + 15552 * (N**3) - 250776 * (N**2) - 239112 * N + 68526

#Px = -81 - 135 * N - 72 * (N**2) - 12 * (N**3)
#Py = 1620 + 1188 * N + 216 * (N ** 2)
#Pz = -3 - Ni

#x = Fraction(Px, Pz)
#y = Fraction(Py, Pz)

#x1 = Fraction(327, 1)
#y1 = Fraction(0, 1)

#print(PointSum(a, b, x, y, x1, y1))
#A = PointSum(a, b, x, y, x1, y1)

x = Fraction(103, 1)
y = Fraction(5824, 1)

#print(Reverse(x, y, Fraction(1), N))

#R = 500
#FindRational(a, b, R)






#data 1
#a = -2
#b = 0
#n = 5
#x = Fraction(-1)
#y = Fraction(1)

#data 2
#a = -5
#b = 5
#n = 5
#x = Fraction(1)
#y = Fraction(1)

n = 30
#print(ScMult(a, b, x, y, n))
#print(GenerateAlotPoints(a, b, x, y))
#print(GenerateNpoints(a, b, x, y, n))
#DrawSumN(a, b, x, y, n)
V = GenerateNpoints(a, b, x, y, n)
xP = V[0]
yP = V[1]
Ans = RevFind(xP, yP, N)
#print(Ans)
a = Ans[0]
b = Ans[1]
c = Ans[2]
ad = a.denominator
bd = b.denominator
cd = c.denominator
a = abs(a.numerator * bd * cd)
b = abs(b.numerator * ad * cd)
c = abs(c.numerator * ad * bd)

print()
print()
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

#print(xP, yP)
#Draw(a, b, xP, yP )
#y = Fraction(-1)
#DrawSum(a, b, [Fraction(1), y], [Fraction(-1),yC(y, a, b)] )
#DrawSum(a, b, [Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)] )
