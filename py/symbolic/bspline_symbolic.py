# import sys
# sys.path.append("..")

import random
import matplotlib.pyplot as plt
import symengine as se

from numpy import array, linspace

from bspline_basis import BsplineBasis

from base_spline import BaseSpline

from deboor import *
t = se.symbols("t")


def latex(x, y):
    return f"$$ x(t) = {se.latex(x)}$$\n\n $$ y(t) = {se.latex(y)}$$\n\n"


def print_latex(x, y):
    print(f"$$ x(t) = {se.latex(x)}$$\n")
    print(f"$$ y(t) = {se.latex(y)}$$")


class BsplineSymbolic(BaseSpline):
    def __init__(self, p: int, points, knots=None, core="rec", closed=False):
        core = "rec"
        """
        Class build p-th degree B-spline at user's points
        spline : C(t) = (x(t), y(t)) where 0 <= t <= 1 
        Construct bspline by basis from Bspline_basis

        Parameters
        ==========
        p : degree of B-spline

        points : array of control points : [[x0, y0], ..., [[xn, yn]]]

        knots : knot vector

        core : type to evalute 
            1. "rec" : Recursive definition B-spline by Cox-De Boor 
            2. "deboor" : Non-recursive de-Boor algorithm 
        """

        super().__init__(p, points, knots, closed)
        self.core = core
        self.spline = ()  # Tuple (x(t), y(t)), x and y are folded piecewise functions
        if core == "deboor":
            self.spline = self.__create_spline_deboor()
        else:
            self.basis = BsplineBasis(self.p, self.knots)
            self.spline = self.create_spline_rec(self.basis, self.points)

        self.func_x, self.func_y = self.get_functions(*self.spline)
        
        
    @staticmethod
    def create_spline_rec(basis, points):
        x = BsplineBasis.zero_piecewise()
        y = BsplineBasis.zero_piecewise()

        for b, p in zip(basis.basis, points):
            x = BsplineBasis.piecewise_sum(1, x, p[0], b)
            y = BsplineBasis.piecewise_sum(1, y, p[1], b)
        return (x, y)

    def __create_spline_deboor(self, func_construct=deBoor):
        args_curve_x = []
        args_curve_y = []
        for k in range(self.p, len(self.knots) - self.p - 1):
            if self.knots[k] != self.knots[k+1]:
                interval = se.Interval(
                    self.knots[k], self.knots[k+1], False, False)
                cond = se.Contains(t, interval)
                curve_part = func_construct(
                    k, t, self.knots, self.points, self.p)
                args_curve_x.append((se.expand(curve_part[0]), cond))
                args_curve_y.append((se.expand(curve_part[1]), cond))

        args_curve_x.append((0, True))
        args_curve_y.append((0, True))
        return (se.Piecewise(*args_curve_x),
                se.Piecewise(*args_curve_y))

    def plot_(self, title="Sym B-spline", N=10000):
        L = self.t_end - self.t_start
        ts = linspace(self.t_start + 1/N, self.t_end - 1/N, N)
        if self.closed:
            plt.plot(self.points[:-self.p][:, 0],
                    self.points[:-self.p][:, 1], 'o--', label='Control Points')
        else:
            plt.plot(self.points[:, 0],
                    self.points[:, 1], 'o--', label='Control Points')
        plt.title(f"B-spline. degree : {self.p}")
        plt.plot(self.func_x(ts)[:-10], self.func_y(ts)[:-10], label=title)
        plt.legend()

    def plot_show(self, title="B-spline", N=10000):
        self.plot_(title, N)
        plt.show()
        plt.close()

    def plot_file(self, path, title="B-spline",  N=10000):
        self.plot_(title, N)
        plt.savefig(path)
        plt.close()

    def get_functions(self, expr_x, expr_y):
        return se.Lambdify(t, expr_x, backend="llvm"), \
            se.Lambdify(t, expr_y, backend="llvm")

    def __call__(self, value):
        return self.func_x(value), self.func_y(value)

    def save_to_latex(self, file="output.tex"):
        with open(file, "w") as f:
            f.write(f"$$ x(t) = {se.latex(self.spline[0])}$$\n")
            f.write(f"$$ y(t) = {se.latex(self.spline[1])}$$")

    def print_latex(self):
        # with * we unpacking tuple spline to spline[0] and spline[1]
        print_latex(*self.spline)


    def create_diff_rec(self):
        p, points, knots = self.construct_data_to_derivative()
        return BsplineSymbolic(p, points, knots, core="rec")

    def create_diff_deboor(self):
        """ Return tuple of Piecewise"""
        p, points, knots = super().construct_data_to_derivative()
        return BsplineSymbolic(p, points, knots, core="deboor")

    def diff(self, method="directly"):
        if self.p == 1:
            raise ValueError("The degree of the spline is already the first!")

        if method == "directly":
            # tuple of Piecewise x'(t), y'(t)
            return se.diff(self.spline[0], t), se.diff(self.spline[1], t)
        else:
            if self.core == "rec":
                return self.create_diff_rec()
            else:
                return self.create_diff_deboor()

def generate_points(n):
    return array([array([i, random.randint(-n,n)]) for i in range(n)])


if __name__ == "__main__":
    """Example creating Bspline using de Boor algorithm"""
    p = 3
    points = [[1, 1],[2, 2], [3, 3],[4, 4]]

    #knots = list(range(len(points) + p + 1))
    spline = BsplineSymbolic(p, points)

    spline.basis.plot_file("basis.png")
    spline.plot_show()
    spline.print_latex()
