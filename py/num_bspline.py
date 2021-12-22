import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate 

import sys
from generate import *

def num_bspline(t, degree, points, knots=None, weights=None):
    i = j = s = l = 0
    n = len(points)
    d = len(points[0])

    assert degree > 1
    assert degree < (n-1)

    if weights == None:
        weights = []
        for i in range(n):
            weights.append(1)

    #generate knot vector
    if knots is None:
        knots = [0]*(degree+1) + [(i)/(n - degree) for i in range(1, n - degree)] + [1]*(degree+1)
    
    assert len(knots) == n + degree + 1
    
    domain = [degree, len(knots) - 1 - degree]

    
    for i in range(domain[0], domain[1]):
        if (t>=knots[i] and t <= knots[i+1]):
            s = i
            break
            
    v = []
    for i in range(n):
        sub_v = []
        for j in range(d):
            sub_v.append(points[i][j] * weights[i])
        sub_v.append(weights[i])
        v.append(sub_v)

    alpha = 0
    for l in range(1, degree + 2):
        for i in range(s, s-degree-1+l, -1):
            den = (knots[i+degree+1-l] - knots[i])
            alpha = (t - knots[i]) / den
            for j in range(d+1):
                v[i][j] = (1 - alpha) * v[i-1][j] + alpha * v[i][j]
    
    res = []
    for i in range(d):
        res.append(v[s][i] / v[s][d])
    return res


def create_curve(p, points, n=10000):
    ts = np.linspace(0, 1, n)
    curve = np.array([num_bspline(t, p, points) for t in ts])
    return curve

def plot_curve(curve, points):
    points = np.array(points)
    curve = np.array(curve)
    plt.scatter(points[:,0], points[:,1], c='g')
    plt.plot(curve[:,0], curve[:,1])



if __name__ == "__main__":    
    p = int(sys.argv[1])  
    points = load_points()
    curve = create_curve(p, points, 1000)
    plot_curve(curve, points)
    plt.show()  
    x, y = curve.T
