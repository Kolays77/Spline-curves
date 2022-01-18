import matplotlib.pyplot as plt
import numpy as np

from generate import *


def length_points(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 -y1)**2)

def compare_curves(xs_1, ys_1, xs_2, ys_2):
    errors = [length_points(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(xs_1, ys_1, xs_2, ys_2)]
    return errors

if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print(f"File1 : {file1}")
    print(f"File2 : {file2}\n")

    points1 = load_points(file1)
    points2 = load_points(file2)

    xs1, ys1 = points1.T
    xs2, ys2 = points2.T


    errors = compare_curves(xs1, ys1, xs2, ys2)
    plt.plot(xs1, ys1, label=file1)
    plt.plot(xs2, ys2, label=file2)
    plt.legend()
    plt.savefig(f"plot_curves.png")
    plt.clf()

    plt.plot(errors)
    plt.yscale("log")
    plt.savefig(f"errors.png")
    plt.clf()