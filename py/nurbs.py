import nurbspy as nrb
import numpy as np
import matplotlib.pyplot as plt

from generate import *

def create_knots(n, p):
    interiors = [(i)/(n - p) for i in range(1, n - p)]
    return np.array([0]*(p+1) + interiors + [1]*(p+1))


def create_curve(p, points, weights, knots, N):
    curve = nrb.NurbsCurve(control_points=points.T, degree=p, knots=knots, weights=weights)
    ts = np.linspace(0,1, 1000)
    points_curve = [curve.get_value(t) for t in ts]
    res = []
    for pp in points_curve:
        res.append([pp[0][0], pp[1][0]])
    return res


def count_dy_dx(p, points, knots, weights, N):
    curve = nrb.NurbsCurve(control_points=points.T, degree=p, weights=weights, knots=knots)
    ts = np.linspace(0,1, 1000)
    der = [curve.get_derivative(t, 1) for t in ts]
    dy_dx = [(d[1]/d[0])[0] for d in der]
    save_vector(dy_dx, f"nurbs_out/dy_dx_num{p}.out")


if __name__ == "__main__":
    points = load_points()
    p = int(sys.argv[1])
    knots = create_knots(len(points), p)
    weights = np.linspace(1.0, 5.0, len(points))
    points_curve = np.array(create_curve(p, points, weights, knots, 1000))
    save_points(points_curve,"nurbs_out/nurbspy_points.out")
    x, y = points_curve.T
    #plt.plot(x,y)
    #plt.show()

    count_dy_dx(p, points, knots, weights, 1000)

