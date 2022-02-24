import sys

import matplotlib.pyplot as plt
import numpy as np

from generate import *


def length_points(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 -y1)**2)

def compare_curves(xs_1, ys_1, xs_2, ys_2):
    errors = [length_points(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(xs_1, ys_1, xs_2, ys_2)]
    return errors

if __name__ == "__main__":

    arg = sys.argv[1]
    file1 = f"out/points_bspline_{arg}.out"
    file2 = f"out/points_num_{arg}.out"

    print(f"File1 : {file1}")
    print(f"File2 : {file2}\n")

    points1 = load_points(file1)
    points2 = load_points(file2)

    xs1, ys1 = points1.T
    xs2, ys2 = points2.T
    cv = load_points("points.in")
    cv_x, cv_y = cv.T
    errors = compare_curves(xs1, ys1, xs2, ys2)
    plt.rcParams["savefig.bbox"] = 'tight'
    plt.rcParams["savefig.dpi"] = 300
    plt.plot(xs1, ys1, label="analytical B-spline")
    plt.plot(xs2, ys2, label="numerical B-spline")
    plt.scatter(cv_x, cv_y, c='g')
    plt.legend()

    plt.savefig(f"src/plot_curves_{arg}.png")
    plt.clf()

    plt.plot(errors)
    plt.yscale("log")
    plt.savefig(f"src/errors_{arg}.png")
    plt.clf()