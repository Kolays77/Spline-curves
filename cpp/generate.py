import numpy as np
from decimal import Decimal
import sys
import math

def generate_points1_array(n):
    return np.array([np.array([i, i]) for i in range(n)])


def generate_points_all_random(n):
    return [[np.random.rand(), np.random.rand()] for i in range(n)]

def generate_points_array_2(n):
    xs = np.sort([np.random.rand() for i in range(n)])
    ys = [np.random.rand() for i in range(n)]
    return np.array( [np.array([x,y]) for x, y in zip(xs, ys)])



def generate_special1(n):
    points1 = generate_points_all_random(n//2)
    points2 = np.array([np.array([i, i]) for i in range(n//2, n)])
    return np.concatenate((points1, points2))

def generate_special2(n):
    points1 = np.array([np.array([i, i]) for i in range(n//2)])
    
    xs = np.sort([n//2 + np.random.rand() for i in range(n)])
    ys = [n* np.random.rand() for i in range(n)]
    
    points2 = np.array( [np.array([x,y]) for x, y in zip(xs, ys)])
    return np.concatenate((points1, points2))


def generate_special1(n):
    points1 = generate_points_all_random(n//2)
    points2 = np.array([np.array([i, i]) for i in range(n//2, n)])
    return np.concatenate((points1, points2))

def generate_points_2(n):
    xs = sorted([np.random.rand() for i in range(n)])
    ys = [np.random.rand() for i in range(n)]
    return [ [x,y] for x, y in zip(xs, ys)]

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

def save_points(points, file_x="cv_x.in", file_y="cv_y.in"):
    file_x = open(file_x, "w")
    file_y = open(file_y, "w")
    for p in points:
        file_x.write(f"{p[0]}\n")
        file_y.write(f"{p[1]}\n")    
    file_x.close()
    file_y.close()
    
if __name__ == "__main__":
    n = int(sys.argv[1])
    points = generate_points_array_2(n)
    save_points(points)