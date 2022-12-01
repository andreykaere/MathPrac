#!/usr/bin/python3

from fractions import Fraction
import matplotlib.pyplot as plt
import math
from math import *
import numpy as np
import itertools
import matplotlib.animation as animation
from time import *

try:
    from point_sum import PointSum
except ModuleNotFoundError:
    from .point_sum import PointSum


def DrawCurve(a, b, rng = 1000):
    x = np.arange(-rng,rng, 0.1)
    yCurveU = []
    yCurveD = []
     
    for i in x.tolist():
        if i**3 + a*i + b >= 1/100:    
            yCurveU += [((i**3 + a * i + b))**(1/2)]
            yCurveD += [-((i**3 + a * i + b))**(1/2)]
        else:
            i0 = int(x.tolist().index(i))
            x[i0] = np.nan
            yCurveU += [np.nan]
            yCurveD += [np.nan]            
   
    yCurveU = np.array(yCurveU)
    yCurveD = np.array(yCurveD)    
    
    plt.plot(x,yCurveU,'b', x, yCurveD, 'b')

# Draws sum of P, Q, coordinates in lists xPoints for P.x, and Q.x 
def DrawSum(a, b, xPoints = [], yPoints = []):

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
    xLine1 = [xr]*100000
    yLine1 = np.arange(-50000, 50000, 1)    
    plt.title('Sum of points')

    plt.plot(xL, yLine, linestyle = '-', linewidth = 1, color = 'y')
    plt.plot(xLine1, yLine1, linestyle = '--', linewidth = 1, color = 'darkmagenta')
    plt.plot(xPoints, yPoints, 'ro', xr, yr, 'go')

# Draws Sum of n Points: P = (x, y), 2*P, 3*P,... n*P
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
            plt.text(x0, y0 + 3, str(i), fontsize=9)
        x0 = ScMult(a, b, x, y, i)[0]
        y0 = ScMult(a, b, x, y, i)[1]
        xP.append(x0)
        yP.append(y0)
    plt.plot([x], [y], 'yo')
    plt.show() 

