import numpy as np
import matplotlib.pyplot as plt
def length_points(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 -y1)**2)

def compare_curves(xs_1, ys_1, xs_2, ys_2, ts=None):
    errors = [length_points(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(xs_1, ys_1, xs_2, ys_2)]        
    return errors

def read_vector(file):
    res = []
    with open(file, "r") as f:
        for line in f:
            res.append(float(line))
    return res

def plot_two_curve(cv_x, cv_y, xs1, ys1, xs2, ys2, label1, label2):
    plt.scatter(cv_x, cv_y, label="points")
    plt.plot(xs1, ys1, label=label1)
    plt.plot(xs2, ys2, label=label2)
    plt.grid()
    plt.legend()
    plt.show()
    
    errors = compare_curves(xs1, ys1, xs2, ys2)
    print(f"Max error {label1} and {label2} : " , max(errors))
    plt.plot(errors, label=f"{label1} vs {label2}")
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    cv_x = read_vector("cv_x.in")
    cv_y = read_vector("cv_y.in")
    xs1 = read_vector("xs_prec.out")
    ys1 = read_vector("ys_prec.out")
    xs2 = read_vector("xs.out")
    ys2 = read_vector("ys.out")
    xs3 = read_vector("num_xs.out")
    ys3 = read_vector("num_ys.out")
    plot_two_curve(cv_x, cv_y, xs2, ys2, xs3, ys3, "bspline.cpp", "num_bspline.cpp")
    plot_two_curve(cv_x, cv_y, xs1, ys1, xs3, ys3, "prec_bspline.cpp", "num_bspline.cpp")
