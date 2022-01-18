import numpy as np
import quadpy
import sys


if __name__ == "__main__":
    n = int(sys.argv[1])
    scheme = quadpy.c1.gauss_legendre(n)
    with open("integrate_weights.txt", "w") as file:
        for w in scheme.weights:
            file.write(f"{w}\n")

    with open("integrate_points.txt", "w") as file:
        for p in scheme.points:
            file.write(f"{p}\n")