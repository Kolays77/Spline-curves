import matplotlib.pyplot as plt
import numpy as np

from bspline import BSpline
from generate import *
from num_bspline import *


def length_points(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 -y1)**2)

def compare_curves(xs_1, ys_1, xs_2, ys_2, ts=None):
    errors = [length_points(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(xs_1, ys_1, xs_2, ys_2)]        
    return errors
   
if __name__ == "__main__":
    file_x1 = sys.argv[1]
    file_y1 = sys.argv[2]
    file_x2 = sys.argv[3]
    file_y2 = sys.argv[4]
    for i in range(4):
        print(sys.argv[i+1])
        
    xs1 = load_vector(file_x1)
    ys1 = load_vector(file_y1)
    xs2 = load_vector(file_x2)
    ys2 = load_vector(file_y2)
    errors = compare_curves(xs1, ys1, xs2, ys2)
    plt.plot(xs1, ys1, label="NURBS")
    plt.plot(xs2, ys2, label="Library NURBS")
    plt.legend()
    plt.show()
    
    plt.plot(errors)
    plt.show()
    
    
    
    
    

    