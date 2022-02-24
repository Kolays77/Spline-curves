import matplotlib.pyplot as plt
import numpy as np

from bspline import BSpline
from generate import *
from num_bspline import *
import sys
     
    
def test_integral(p, n=None):
    points = load_points()
    spl = BSpline(p, points, closed=False)
    spl.plot(fill=True)
    plt.show()
    N = 1000
    integral1 = spl.num_integral(1/N)
    integral2 = spl.integral(N)
    
    print(f"num {integral1=}")
    print(f"formulae {integral2=}")
 

def test_derivative(p, n=None):
    points = load_points("points.in")
    spl = BSpline(p, points, closed=False)
    spl.plot()
    plt.show()
    N = 1000
    with open("dy_dx.out", "w") as out:
        for t in np.linspace(spl.t_start, spl.t_end, N):
            out.write(f"{spl.dy_dx(t)}\n")

    der_spl = spl.create_diff_bspline()
    der_spl.print_coef(file_x="der_x_coef.out", 
                        file_y="der_y_coef.out")
    
    dy_dx = []
    for k, coef in zip(der_spl.intervals, der_spl.curve_coefs):
        ts = np.linspace(spl.knots[k], spl.knots[k+1], N)
        temp = [coef[1](t)/coef[0](t) for t in ts]
        dy_dx.append(temp)
        
p = int(sys.argv[1])
n = int(sys.argv[2])

test_integral(p, n)