import sys
from itertools import accumulate

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import scale
from numpy import array, concatenate, linspace, random
from scipy import integrate

from base_spline import BaseSpline
from deboor import *
from generate import *


def length_points(P0, P1):
    return np.sqrt((P0[0] - P1[0])**2 + (P0[1] - P1[1])**2)


def check_first_last_point():
    print(f"First point")


class BSpline(BaseSpline):
    def __init__(self, p: int, points, knots=None, weights=None, closed=False):
        super().__init__(p, points, knots, weights, closed)
        self.intervals = self.create_intervals()
        self.coefs = []

        if self.dim > 1:
            for k in self.intervals:
                self.coefs.append(deBoor_coef(
                    k, self.knots, self.points, self.p))
        else:
            for k in self.intervals:
                self.coefs.append(deBoor_coef_1d(
                    k, self.knots, self.points, self.p))
        if self.dim > 1:
            _, self.der_coefs, _ = self._derivative_directly()


    def __call__(self, arg):
        for k, coef in zip(self.intervals, self.coefs):
            if self.knots[k] <= arg <= self.knots[k+1]:
                return [coef[i](arg) for i in range(self.dim)]

    def create_intervals(self):
        intervals = []
        for k in range(self.p, len(self.knots) - self.p - 1):
            if self.knots[k] != self.knots[k+1]:
                intervals.append(k)
        return intervals

    def get_intervals(self):
        res_interval_pairs = []
        for i in self.intervals:
            res_interval_pairs.append([self.knots[i],self.knots[i+1]])
        return res_interval_pairs
    
    def get_curve_points(self, N):
        t0 = self.t_start
        tn = self.t_end
        ts = linspace(t0, tn, N)      
        points = array([self.__call__(t) for t in ts])     
        
        if self.dim > 1:         
            return points[:, 0], points[:, 1]
        else:
            return points
        
    def plot(self, N=1000, fill=False):
        xs, ys = self.get_curve_points(N)
        plt.scatter(array(self.points)[:, 0], array(self.points)[:, 1])
        plt.plot(xs, ys, label="B-spline")
        if fill:
            plt.fill_between(xs, ys, 0)
        plt.legend()


    def print_coef(self, der=False):
        if self.dim > 1:    
            if der:
                coefs = self.der_coefs
                x_out = open("der_x_coef.out", "w")
                y_out = open("der_y_coef.out", "w")
            else:
                coefs = self.coefs   
                x_out = open("x_coef.out", "w")
                y_out = open("y_coef.out", "w")
                
            for k, coef in zip(self.intervals, coefs):
                x_out.write(
                    f"[{self.knots[k]}, {self.knots[k+1]}] {coef[0].coef}\n")
                y_out.write(
                    f"[{self.knots[k]}, {self.knots[k+1]}] {coef[1].coef}\n")        
            x_out.close()
            y_out.close()
        else:
            out = open("coefs.out", "w")
            for k, coef in zip(self.intervals, self.coefs):
                out.write(f"[{self.knots[k]}, {self.knots[k+1]}] {coef[0].coef}\n")
            out.close()


    ### DERIVATIVES     
    
    def create_diff_bspline(self):
        p, points, knots = super().construct_data_to_derivative()
        return BSpline(p, points, knots)
    
    def _derivative_directly(self):
        coef_der = []
        for coef in self.coefs:
            coef_der.append([coef[0].deriv(), coef[1].deriv()])
        return (self.p-1, coef_der, self.get_intervals())

    def dy_dx(self, arg, N=1000):
        # dy/dx = y'(t)/x'(t)
        for k, coef, der_coef in zip(self.intervals, self.coefs, self.der_coefs):
            if self.knots[k] <= arg <= self.knots[k+1]:
                x, y = der_coef[0](arg), der_coef[1](arg)
                return y/x
    
    ### INTEGRALS
    
    def integral(self, N=1000, t0=0, t1=1):
        # integral :  integral{t0}{t1}{y(t)*x'(t)dt}    
        res_sum = 0   
        for k, coef, der_coef in zip(self.intervals, self.coefs, self.der_coefs):
            F = np.poly1d(der_coef[0]*coef[1]).integ() # Antiderivative poly
            if self.knots[k+1] > t0 and self.knots[k] < t1:
                start = max(t0, self.knots[k])
                finish = min(t1, self.knots[k+1])
                res_sum += F(finish) - F(start)            
        return res_sum
 
    def num_integral(self, h):
        xs, ys = self.get_curve_points(int(1/h))
        return integrate.simpson(ys, xs)
    
        
if __name__ == "__main__":
    if len(sys.argv) > 1:    
        p = int(sys.argv[1])
        points = load_points()     
        curve = BSpline(p, points, knots=None, weights=None)
        curve.plot()
        xs, ys = curve.get_curve_points(1000)
        save_vector(xs, "xs.out")
        save_vector(ys, "ys.out")

        plt.show()
        