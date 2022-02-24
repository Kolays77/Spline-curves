from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

def deBoor(k, t, knots, points, p):
    """
    Evaluates S(x).
    Parameters
    ==========
    k: index of knot interval that contains t
    t: position. OR Symbol "t" then the output will be the expression
    knots: array of knot positions, needs to be padded as described above
    points: array of control points
    p: degree of B-spline
    """
    d = [points[j + k - p] for j in range(0, p+1)]
    for r in range(1, p+1):
        for j in range(p, r-1, -1):
            alpha = (t - knots[j+k-p]) / (knots[j+1+k-r] - knots[j+k-p])
            d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j]
    return d[p]


def deBoor_coef(k, knots, cv, p):
    #With optimization
    # coefs : c_p*t^p + ... + c_0
    coefs =[ [np.poly1d(cv[j + k - p][0]), 
               np.poly1d(cv[j + k - p][1])] for j in range(0, p+1)]
    dim = len(cv [0]) # (x, y) -> 2, (x, y, z ) -> 3
    for r in range(1, p+1):
        for j in range(p, r-1, -1):
            den = (knots[j+1+k-r] - knots[j+k-p])
            A = 1 / den
            B = - knots[j+k-p] / den
            C = -A
            D = knots[j+1+k-r] / den
            
            for d in range(dim):       
                temp1 = A * coefs[j][d] * np.poly1d([1, 0]) + B * coefs[j][d]
                temp2 = C * coefs[j-1][d] * np.poly1d([1, 0]) + D * coefs[j-1][d]
                coefs[j][d] = np.poly1d(temp1 + temp2) 
    return [coefs[p][0], coefs[p][1]]

def deBoor_coef_1d(k, knots, cv, p):
    #With optimization
    # coefs : c_p*t^p + ... + c_0
    coefs =[ np.poly1d(cv[j + k - p]) for j in range(0, p+1)]
    for r in range(1, p+1):
        for j in range(p, r-1, -1):
            den = (knots[j+1+k-r] - knots[j+k-p])
            A = 1 / den
            B = - knots[j+k-p] / den
            C = -A
            D = knots[j+1+k-r] / den
            temp1 = A * coefs[j] * np.poly1d([1, 0]) + B * coefs[j]
            temp2 = C * coefs[j-1] * np.poly1d([1, 0]) + D * coefs[j-1]
            coefs[j] = np.poly1d(temp1 + temp2) 
    
    return [coefs[p]]