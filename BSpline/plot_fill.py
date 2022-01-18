import matplotlib.pyplot as plt
import numpy as np

from generate import *


if __name__ == "__main__":
    file = sys.argv[1]
    print(f"File : {file}")

    points = load_points(file)
    cv = load_points("points.in")

    xs, ys = points.T
    cv_x, cv_y = cv.T

    plt.scatter(cv_x, cv_y, c='g')
    plt.plot(xs, ys, label=file, c='r')
    plt.fill_between(xs, np.zeros(len(points)), ys)
    plt.legend()
    plt.show()

    #plt.savefig(f"plot_fill.png")
