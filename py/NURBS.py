from os import fdopen
import sys
sys.path.append('..')

from base_spline import BaseSpline
from bspline import BSpline
import matplotlib.pyplot as plt
import numpy as np
from generate import *
from compare import *
 
class NURBS(BaseSpline):
    def __init__(self, p: int, points, knots=None, weights=None, closed=False):
        super().__init__(p, points, knots, weights, closed)
        self.num = BSpline(self.p, points, knots, self.weights)
        self.den = BSpline(self.p, np.ones(self.n), knots, self.weights)
        self.print_coef_num()
        self.print_coef_den()
        
    def print_coef_den(self):
        den_coef = open("den_coef.out", "w")
        for k, den in zip(self.den.intervals, self.den.coefs):
            den_coef.write(f"[{self.knots[k]}, {self.knots[k+1]}] {den[0].coef}\n")
        den_coef.close()
        
    def print_coef_num(self):
        num_coef_x = open("num_coef_x.out", "w")
        num_coef_y = open("num_coef_y.out", "w")
        
        for k, num in zip(self.num.intervals, self.num.coefs):
            num_coef_x.write(f"[{self.knots[k]}, {self.knots[k+1]}] {num[0].coef}\n")
            num_coef_y.write(f"[{self.knots[k]}, {self.knots[k+1]}] {num[1].coef}\n")
            
        num_coef_y.close()
        num_coef_x.close()
        
    
    def get_curve_points(self, N=1000):    
        points = []
        i = 0
        ts = np.linspace(self.t_start, self.t_end, N)
        t = ts[0]
        for k, num, den in zip(self.num.intervals, self.num.coefs, self.den.coefs):
            while(t < self.knots[k+1]):
                temp  = [num[0](t) / den[0](t), 
                         num[1](t) / den[0](t)]
                points.append(temp)
                i += 1
                t = ts[i]
        return np.array(points)       
        
    def plot(self, N=1000):
        xs, ys = self.get_curve_points(N).T
        #plt.scatter(self.points[:,0], self.points[:,1])
        plt.plot(xs, ys, label="NURBS")
        plt.legend()
    
    
def test_circle():
    p = 2
    points = np.array([[1,0], [1,1],[0,1],[-1,1],[-1,0], [-1,-1], [0,-1], [1,-1], [1,0]])
    knots = np.array([0]*(p+1) +[np.pi/2, np.pi/2, np.pi, np.pi, 1.5*np.pi, 1.5*np.pi] + [2*np.pi]*(p+1))
    weights = [1, np.sqrt(2)/2]*4 + [1]
    curve = NURBS(p, points, knots=knots, weights=weights)
    N = 1000
    ts = np.linspace(0, 2*np.pi, 4*N)
    plt.plot(np.cos(ts), np.sin(ts), label="Circle")
    cv_x, cv_y = points.T
    plt.scatter(cv_x, cv_y)
    curve.plot()
    plt.legend()
    plt.show() 
    
    x1, y1 = curve.get_curve_points(N).T
    x2, y2 = np.cos(ts), np.sin(ts)
    
    
    save_vector(y1, "y1.out")
    save_vector(y2, "y2.out")
    
    errors = compare_curves(x1, y1, x2, y2)
    plt.plot(errors)
    plt.show()
   

def test():
    p = int(sys.argv[1])
    points = load_points()   
    w = list(range(1, len(points)+1))
    w = [7, 13, 7, 3, 2, 1, 8, 12, 11, 9]
    curve = NURBS(p, points, knots=None, weights=w)
    curve.plot()
    
    cv_x, cv_y = points.T
    plt.scatter(cv_x, cv_y)
    plt.show()

if __name__ == "__main__":
    test()
    