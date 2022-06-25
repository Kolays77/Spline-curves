import numpy as np
import sys


def generate_square_for_integral(n, size):
    xs = np.sort([size*np.random.rand() for i in range(n)])
    ys = [size*np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])


def generate_points1_array(n):
    return np.array([np.array([i, np.random.randint(-n,n)]) for i in range(n)])

def generate_points_all_random(n):
    xs = [np.random.rand() for i in range(n)]
    ys = [n*np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])

def generate_points(n):
    xs = np.sort([n*np.random.rand() for i in range(n)])
    ys = [n/2.0 - n*np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])

def generate_points3(n):
    return np.array([np.array([i, i]) for i in range(n)])

def load_points(file="points.in"):
    points = []
    with open(file, "r") as file:
        for line in file:
            point = np.array(list(map(float, line.split())))
            points.append(point)
    return np.array(points)

def save_points(points, file="points.in"):
    with open(file, "w") as file:
        for p in points:
            file.write(f"{p[0]}\t{p[1]}\n")

def save_vector(vec, file):
    with open(file, "w") as file:
        for p in vec:
            file.write(f"{p}\n")

def load_vector(file):
    vec = []
    with open(file, "r") as file:
        for line in file:
            vec.append(float(line ))
    return np.array(vec)


if __name__ == "__main__":
    n = int(sys.argv[1])
    points = generate_points(n)
    save_points(points, "points.in")
    
