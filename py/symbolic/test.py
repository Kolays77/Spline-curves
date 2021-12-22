from bspline.py.symbolic.bspline_symbolic import BsplineSymbolic
import os
import shutil
import numpy as np
import sys
import matplotlib.pyplot as plt
import random 
import symengine as se
from os import path
cores = ["rec", "deboor"]
SYSTEM = "UNIX"

#SYSTEM = "WINDOWS" # uncomment this for Windows system

if SYSTEM == "WINDOWS":
    os.sep = '\\'
else:
    os.sep = '/'
    
#TODO do it for float and more complicated 
def generate_points(n):
    return np.array([np.array([i, random.randint(-n,n)]) for i in range(n)])

def length_points(P0, P1):
    return np.sqrt( (P0[0] - P1[0])**2 + (P0[1] - P1[1])**2)
    

def create_knots1(p, n):
    interiors = [se.Rational(i,(n - p)) for i in range(1, n - p)]
    return [0]*(p+1) + interiors + [1]*(p+1)

def create_knots2(p, n):
    interiors = [se.Rational(1, 8) for _ in range(1, n - p)]
    return [0]*(p+1) + interiors + [1]*(p+1)

def get_points():
    return np.array([[0,0], [1,5], [2,-5], [3,0], [4,-6],[5,-2]])

import random
def test_curve_from_knot_vector():
    p = 3
    n = 6
    m = p + n + 1
    points = get_points()
    
    knot_uniform = [i for i in range(m)]
    knot_opened1 = create_knots1(p, n)
    knot_opened2 = create_knots2(p, n)
    knot_random = sorted([random.randint(0, m) for _ in range(m)])
    dir_path =  path.join("Tests", "test_curve_from_knot_vector")
    
    for knots, name in zip([knot_uniform, knot_opened1, knot_opened2, knot_random], 
                           ["uniform", "opened1", "opened2", "knot_random"]):
        
        spline = BsplineSymbolic(p, points, knots)
        plot_basis_path = path.join(dir_path, name + "_basis")
        plot_curve_path = path.join(dir_path, name)
        spline.basis.plot_file(plot_basis_path, str(knots), 10000)
        spline.plot_file(plot_curve_path, N=10000)

def test_open_closed():
    dir_path = path.join("Tests","test_open_closed")
    
    for p in range(3, 10):
        points = generate_points(p + 3)
        spline_closed = BsplineSymbolic(p, points, closed=True)    
        spline_opened = BsplineSymbolic(p, points, closed=False)    
        
        opened_path = path.join(dir_path, f"_opened_{p}")
        closed_path = path.join(dir_path, f"_closed_{p}")
        
        spline_closed.plot_file(opened_path)
        spline_opened.plot_file(closed_path)
        


def test_spline_degree():
    points = generate_points(11)
    degrees = [2,3,4,5,6,7,8,9]
    dir_path =  path.join("Tests", "test_spline_degree")
    for p in degrees:
        spline = BsplineSymbolic(p, points, core="rec")
        file_path = path.join(dir_path, f"degree_{p}.tex")
        file = open(file_path, "w+")
        file.write(f"SPLINE : \n\n")
        file.write(se.latex(*spline.spline))
        spline.plot_file(path.join(dir_path, f"degree_{p}"))


        
def test_compare_splines():
    """
    Comparing my recursive Bspline by definition
    and optimized non-recursive symbolic deBoor algorithm
    """
    ERROR = 1e-8
    max_error = 0
    p = 5
    points = generate_points(10)
    spline_rec = BsplineSymbolic(p, points, core="rec", closed=False)
    spline_deboor = BsplineSymbolic(p, points, core="deboor", closed=False)
        
    ts = np.linspace(0, 1, 10000)
    max_error = np.max([length_points(spline_rec(t), spline_deboor(t))  for t in ts])
    
    dir_path =  path.join("Tests", "test_compare_splines")
    file_path = path.join(dir_path,  f"Compare_splines.tex")
    file = open(file_path, "w+")
    
    file.write(f"SPLINE method Recursive \n\n")    
    file.write(se.latex(*spline_rec.spline))
    
    file.write(f"SPLINE method de Boor \n\n")    
    file.write(se.latex(*spline_deboor.spline))
    
    file.write(f"Max error by methods = {max_error}\n")
    
    file.close()
    
    plt.scatter(spline_rec.points[:, 0], spline_rec.points[:, 1])
    for spline in [spline_rec, spline_deboor]:    
        plt.plot(spline.points[:, 0], spline.points[:, 1])
        plt.plot(spline.func_x(ts), spline.func_y(ts))    
        plt.savefig(path.join(dir_path, "Recursive and de Boor"))
    plt.close()
    
    # Error accumulation
    if max_error >= ERROR:
        raise ValueError(f"Max error {max_error} more then {ERROR}.")

    

def test_derivatives():
    p = 6
    points = generate_points(10)    
    for core in cores:
        spline = BsplineSymbolic(p, points, core=core)
        dir_diffs = spline.diff()
        spline_der = spline.diff(core)
        path_file = path.join("Tests", "test_derivatives", f"Test_derivative_method_{core}.tex")
        file = open(path_file, "w+")
        file.write(f"SPLINE method {core} : \n\n")    
        file.write(se.latex(*spline.spline))
        file.write(f"DIRECTLY derivative : \n\n")
        file.write(se.latex(*dir_diffs))
        file.write(f"Formula derivative : \n\n")
        file.write(se.latex(*spline_der.spline))
        file.close()

tests = [test_compare_splines, 
             test_derivatives, 
             test_spline_degree,
             test_curve_from_knot_vector, 
             test_open_closed]
    
    
def mkdir_test(name):
    os.mkdir(f"Tests/{name}")

def delete_tests():
    try:
        shutil.rmtree("Tests")
    except Exception:
        pass

def run_tests():
    delete_tests()
    os.mkdir("Tests")
    for test in tests:    
        try:
            mkdir_test(test.__name__)
            test()
        except Exception as err:
            print(f"{test.__name__}  : NOT OK : {err}")
        else:
            print(f"{test.__name__} : OK") 
    
if __name__ == "__main__":
    """
    Usage : 
        1. python3 test.py run - Running all tests with folder creation.
        2. python3 test.py del - Deleting all test folders.
    """
    #run_tests()
    if len(sys.argv) == 1:
        print("You must use the command line arguments. See Usage test.py")
    else:
        if sys.argv[1] == "run":
            run_tests()
        if sys.argv[1] == "delete":
            delete_tests() 