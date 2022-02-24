import numpy as np
import nurbspy as nrb
import matplotlib.pyplot as plt
import sys
from generate import *

def create_knots(n, p):
    interiors = [(i)/(n - p) for i in range(1, n - p)]
    return np.array([0]*(p+1) + interiors + [1]*(p+1))

if __name__ == "__main__":
    points = load_points()
    p = int(sys.argv[1])
    knots = create_knots(len(points), p)
    curve = nrb.NurbsCurve(control_points=points.T, degree=p, knots=knots)
    
    ts = np.linspace(0,1, 1000)
    points_curve = np.array([curve.get_value(t) for t in ts])
    
    xs = []
    ys = []
    for p in points_curve:    
        xs.append(p[0][0])
        ys.append(p[1][0])
    
    save_vector(xs, "xs_nurbspy.out")
    save_vector(ys, "ys_nurbspy.out")

    ders = np.array([ curve.get_derivative(t, order=1) for t in ts])
    der_x = []
    der_y = []
    for d in ders:    
        der_x.append(d[0][0])
        der_y.append(d[1][0])
    
    dy_dx = [ y / x for x, y in zip(der_x, der_y)]
    
    save_vector(der_x, "dy_dx_nurbspy.out")
    save_vector(der_x, "der_x_nurbspy.out")
    save_vector(der_y, "der_y_nurbspy.out")
    
    plt.scatter(points[:,0], points[:,1], label="CV")
    plt.plot(xs, ys, label="NURBSPY")
    plt.legend()
    plt.show()