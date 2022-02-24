import sys

import matplotlib.pyplot as plt
import numpy as np

from generate import *


if __name__ == "__main__":
    file_cv = "points.in"
    file = sys.argv[1]
    if len(sys.argv) > 2:
        file_cv = sys.argv[2]

    points = load_points(file)
    cv = load_points(file_cv)

    xs, ys = points.T
    cv_x, cv_y = cv.T

    plt.scatter(cv_x, cv_y, c='g')
    plt.plot(xs, ys, label=file, c='b')
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["savefig.bbox"] = 'tight'

    if len(sys.argv) > 2:
        plt.fill_between(xs, np.zeros(len(points)), ys)

    plt.legend()
    plt.show()

    plt.savefig(f"plot.png")
