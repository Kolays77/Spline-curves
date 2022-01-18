import nurbspy as nrb

from generate import *

def create_knots(n, p):
    interiors = [(i)/(n - p) for i in range(1, n - p)]
    return np.array([0]*(p+1) + interiors + [1]*(p+1))


def create_curve(p, points, weights, knots, N):
    curve = nrb.NurbsCurve(control_points=points.T, degree=p, knots=knots)
    ts = np.linspace(0,1, 1000)
    points_curve = [curve.get_value(t) for t in ts]
    res = []
    for pp in points_curve:
        res.append([pp[0][0], pp[1][0]])
    return res

if __name__ == "__main__":
    points = load_points()
    weights = np.linspace(1.0, 5.0, len(points))
    p = int(sys.argv[1])
    knots = create_knots(len(points), p)
    points_curve = create_curve(p, points, weights, knots, 1000)
    save_points(points_curve,"nurbspy_points.out")

