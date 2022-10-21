#!/usr/bin/python3

from fractions import Fraction
import matplotlib.pyplot as plt
import math
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
        m = (3 * xp * xp + a) / (2 * yp)
        xr = m * m - xp - xq
        yr = yp + m * (xr - xp)           

    else:
        m = (yp - yq)/ (xp - xq)
        xr = m * m - xp - xq
        yr = yp + m * (xr - xp)   
    
    return [xr, -yr]

def ScMult(a, b, x, y, n):
    x0 = x
    y0 = y
    for i in range(n - 1):
        A = PointSum(a, b, x, y, x0, y0)
        x = A[0]
        y = A[1]
    return [x, y]

def Draw(a = -5, b = 5, xPoints = [], yPoints = []):
    x = np.arange(-10,10,0.01)
    yu = ((x**3 + a * x + b))**(1/2)
    yd = - ((x**3 + a * x + b))**(1/2)
    plt.title('My first plot')
    plt.plot(x,yu,'b', x, yd, 'g', xPoints, yPoints, 'ro')
    plt.show()    

# to draw plt.show()
def DrawSum(a = -5, b = 5, xPoints = [], yPoints = []):
    x = np.arange(-10,10,0.01)
    xL = np.arange(-5,5,0.1)
    # vertex = b**(1/2)
 
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
    
    #yCurveU = ((x**3 + a * x + b))**(1/2)
    #yCurveD = - ((x**3 + a * x + b))**(1/2)
   
    xp = xPoints[0]
    yp = yPoints[0]
    xq = xPoints[1]
    yq = yPoints[1] 
    A = PointSum(a, b, xp, yp, xq, yq)
    #print(A)
    xr = A[0]
    yr = A[1]
    xr = [xr]
    yr = [yr]
    
    if (xp == xq) and (yp == yq):
        m = (3 * xp * xp + a) / (2 * yp)
    
    else:
        m = (yp - yq)/ (xp - xq)
    yLine = m * (xL - xq) + yq
    xLine1 = [xr]*2000
    yLine1 = np.arange(-10, 10, 0.01)    
    plt.title('Sum of two points')

    plt.plot(x,yCurveU,'b', x, yCurveD, 'b')
    plt.plot(xL, yLine, linestyle = '-', linewidth = 1, color = 'y')
    plt.plot(xLine1, yLine1, linestyle = '--', linewidth = 1, color = 'darkmagenta')
    plt.plot(xPoints, yPoints, 'ro', xr, yr, 'go')
    #plt.show()

def yC(x, a, b):
    y2 = (x**3 + a * x + b)
    #y2 = -y2
    y2n =  abs(y2.numerator)
    y2d = abs(y2.denominator)
    if (abs(y2n) > 10**10) or (abs(y2d) > 10**10):
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
    fig, ax = plt.subplots(figsize=(9, 9))
    
    
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

#data 1
#a = -2
#b = 0
#n = 5
#x = Fraction(-1)
#y = Fraction(1)

#data 2
a = -5
b = 5
n = 5
x = Fraction(1)
y = Fraction(1)


print(ScMult(a, b, x, y, n))
DrawSumN(a, b, x, y, n)
#print(xP, yP)
#Draw(a, b, xP, yP )
#y = Fraction(-1)
#DrawSum(a, b, [Fraction(1), y], [Fraction(-1),yC(y, a, b)] )
#DrawSum(a, b, [Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)] )
