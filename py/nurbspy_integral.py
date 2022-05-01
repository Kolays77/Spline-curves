import numpy as np
import nurbspy as nrb
import matplotlib.pyplot as plt
import sys

import scipy.integrate as integrate

def create_knots(n, p):
    interiors = [(i)/(n - p) for i in range(1, n - p)]
    return np.array([0.0]*(p+1) + interiors + [1.0]*(p+1))


def generate_square_for_integral(n, size):
    xs = np.sort([size*np.random.rand() for i in range(n)])
    ys = [size*np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])

def load_points(file="points.in"):
    points = []

    with open(file, "r") as file:
        points = np.array([np.array(list(map(float, line.split()))) for line in file])    
    return np.array(points)



if __name__ == "__main__":
    p = int(sys.argv[1])
    n = int(sys.argv[2])

    points = load_points("points.txt")
    knots = create_knots(n, p)
    w = np.linspace(1.0, 2.0, n)

    curve = nrb.NurbsCurve(control_points=points.T,  weights=w, degree=p, knots=knots)
    
    ts = np.linspace(0.0, 1.0, 1000)
    points_curve = np.array([curve.get_value(t) for t in ts])

    xs = []
    ys = []
    for p in points_curve:    
        xs.append(p[0][0])
        ys.append(p[1][0])

    ders = np.array([ curve.get_derivative(t, order=1) for t in ts])
    der_x = []
    der_y = []
    for d in ders:    
        der_x.append(d[0][0])
        der_y.append(d[1][0])

    
    ys = np.array(ys)
    der_x = np.array(der_x)
    
    f = np.array(ys * der_x)
    int_simp = integrate.simpson(f, ts)
    print("INTEGRAL : " , int_simp)
    
    # plt.scatter(points[:,0], points[:,1], label="CV")
    # plt.plot(xs, ys, label="NURBSPY")
    # plt.legend()
    # plt.show()