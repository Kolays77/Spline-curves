import sys

import matplotlib.pyplot as plt
from generate import *


if __name__ == "__main__":
    points = load_points(sys.argv[1])
    xs, ys = points.T
    xs = [int(x) for x in xs]
    plt.plot(xs, ys)
    plt.grid()
    plt.rcParams["savefig.bbox"] = 'tight'
    plt.rcParams["savefig.dpi"] = 300
    plt.yscale("log")
    plt.savefig("errors")
