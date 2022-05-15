import matplotlib.pyplot as plt
import sys
import numpy as np


plt.rcParams["savefig.bbox"] = 'tight'
plt.rcParams["savefig.dpi"] = 300

errors = []

path = sys.argv[1]  
with open(path) as file:
    for line in file:
        errors.append(float(line))

x = np.arange(0, len(errors), 1)

plt.scatter(x, errors);
plt.yscale("log")
plt.show()
