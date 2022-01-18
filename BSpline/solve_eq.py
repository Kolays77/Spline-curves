import scipy
import numpy as np

if __name__ == "__main__":
    poly = np.poly1d([1,2,3,4,5])
    np.set_printoptions(precision=15)
    print(poly.r)