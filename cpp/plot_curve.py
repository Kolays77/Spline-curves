import matplotlib.pyplot as plt
import sys
def get_vec_from_file(file):
    vec = []
    with open(file, "r") as out:
        for v in out:
            vec.append(float(v))
    return vec


if __name__ == "__main__":
    xs = get_vec_from_file(sys.argv[1])
    ys = get_vec_from_file(sys.argv[2])
    cv_x = get_vec_from_file("cv_x.in")
    cv_y = get_vec_from_file("cv_y.in")
    plt.plot(cv_x, cv_y, "go--")
    plt.plot(xs, ys)
    plt.show()