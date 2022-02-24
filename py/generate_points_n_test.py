import sys

from generate import *

n_test = int(sys.argv[1])
n_points = int(sys.argv[2])
size = float(sys.argv[3])
dir = sys.argv[4]

for i in range(1, n_test+1):
    points = generate_square_for_integral(n_points, size)
    file = dir+"points_"+str(i)+".in"
    save_points(points, file)
