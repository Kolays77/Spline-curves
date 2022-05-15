import functools
import symengine as se
import matplotlib.pyplot as plt
from numpy import array, linspace
import matplotlib as mpl

t = se.symbols("t") 

class BsplineBasis:
    def __init__(self, p: int, knots: list):
        """
        Class build set basis functions of B-spline 
        """
        self.p = p
        self.knots = array(knots)
        self.m = len(knots) - p - 1

        self.basis = self.create_basis()

    def _divide(self, num, den):
        if den == 0:
            return 0
        else:
            return num / den

    @staticmethod
    def zero_piecewise():
        return se.Piecewise((0, True))


    @staticmethod
    def piecewise_sum(a, pw1, b, pw2):
        """ construct  pw : a * pw1 + b*pw2"""
        dict1, dict2 = {}, {}
        # Last args are zero_piecewise, skip this.
        # Add zero_piecewise at the end of the function
        for i in range(0, len(pw1.args)-2, 2):
            dict1[pw1.args[i+1]] = pw1.args[i]*a
        for i in range(0, len(pw2.args)-2, 2):
            dict2[pw2.args[i+1]] = pw2.args[i]*b
        lst = []
        for key in dict2:
            if key in dict1:
                dict1[key] += dict2[key]
            else:
                dict1[key] = dict2[key]
        for key in dict1:
            lst.append((se.expand(dict1[key]), key))
        # Always at the end zero_piecewise
        lst.append((0, True))
        return se.Piecewise(*lst)

    @functools.cache
    def _B(self, i: int, p: int):
        if p == 0:
            if self.knots[i] == self.knots[i+1]:
                return self.zero_piecewise()
            else:
                return se.Piecewise((1, se.Contains(t, se.Interval(self.knots[i], self.knots[i+1], False, False))),
                                    (0, True))

        sum = self.piecewise_sum(self._divide(t-self.knots[i], self.knots[i+p] - self.knots[i]),
                            self._B(i, p-1),
                            self._divide(
                                self.knots[i+p+1]-t, self.knots[i+p+1] - self.knots[i+1]),
                            self._B(i+1, p-1))
        return sum

    def create_basis(self):
        basis = []
        for i in range(self.m):
            b = self._B(i, self.p)
            basis.append(b)
        return basis

    def plot_(self, title=None, N=1000):
        ts = linspace(self.knots[0], self.knots[-1], N)
        for i, b in enumerate(self.basis):
            plt.plot(ts, se.Lambdify(t, b, backend="llvm")(ts), label=f"$B_" + "{" f"{i}p" + "}$")
        if title != None:
            plt.title(title)
        plt.legend()
    
    def plot_file(self, path="Basis", title=None,N=1000):
        self.plot_(title, N)
        plt.savefig(path)
        plt.close()

    def plot_show(self,title=None, N=1000):
        self.plot_(title,N)
        plt.show()


if __name__ == "__main__":
    # n = 4, p = 3, m = 8
    mpl.rc('text', usetex = True)

    knots = [0, 1, 2, 3, 4, 5, 6, 7]    
    basis_obj = BsplineBasis(3, knots)
    basis_obj.plot_file("basis.png")
