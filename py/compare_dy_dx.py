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

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    arg = sys.argv[3]

    print(f"File1 : {file1}")
    print(f"File2 : {file2}\n")

    vec1 = load_vector(file1)
    vec2 = load_vector(file2)
    plt.rcParams["savefig.bbox"] = 'tight'
    plt.rcParams["savefig.dpi"] = 300
    errors = [ np.abs(v1-v2) for v1, v2 in zip(vec1, vec2)]
    plt.plot(errors)
    plt.yscale("log")
    plt.savefig(f"nurbs_out/dy_dx_errors_{arg}.png")
    plt.clf()