import numpy as np
from decimal import Decimal
import sys
import math

def generate_points1_array(n):
    return np.array([np.array([i, np.random.randint(-n,n)]) for i in range(n)])

def generate_points1(n):
    return [[i, np.random.randint(-n,n)] for i in range(n)]

def generate_points(n):
    xs = np.sort([np.random.rand() for i in range(n)])
    ys = [n*np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])

def generate_points3(n):
    return np.array([np.array([i, i]) for i in range(n)])


def generate_full_circle(n):
    return np.array([np.array([math.cos(np.pi*2*i/n), math.sin(np.pi*2*i/n)]) for i in range(n)])


def generate_half_circle(n):
    return np.array([np.array([np.cos(i*np.pi/(n-1)),  
                               np.sin(i*np.pi/(n-1))]) for i in range(0, n)])
    
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
    points = generate_points1_array(n)
    save_points(points, "points.in")