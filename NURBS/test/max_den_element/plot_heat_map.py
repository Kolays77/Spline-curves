
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
import pandas as pd
import sys
from matplotlib.colors import LogNorm, Normalize

# usage python3 plot_heat_map.py matrix_1_2_10_12_100_random.out matrix_2_2_10_12_100_random.out 

plt.rcParams["savefig.bbox"] = 'tight'
plt.rcParams["savefig.dpi"] = 300

p_start = 0
p_end = 0
n_start = 0
n_end = 0

val_1 = []
with open(sys.argv[1]) as file:
    params = list(map(int, file.readline().split()))
    p_start = params[0]
    p_end = params[1]
    n_start = params[2]
    n_end = params[3]
    for line in file:
        l = list(map(float, line.split()))
        val_1.append(l)
val_1 = np.array(val_1)

val_2 = []
with open(sys.argv[2]) as file:
    params = list(map(int, file.readline().split()))
    for line in file:
        l = list(map(float, line.split()))
        val_2.append(l)

ps = list(range(p_start, p_end + 1))
ns = list(range(n_start, n_end + 1))
val_2 = np.array(val_2)

data1 = pd.DataFrame(val_1, columns=ns, index=ps)
ax = sns.heatmap(data1, norm=LogNorm(), cmap= 'coolwarm', vmin = val_1.min(), vmax = val_1.max())
plt.xlabel("Количество точек")
plt.ylabel("Степень кривой")
plt.savefig(f"plot_1_{p_start}_{p_end}_{n_start}_{n_end}")
plt.close()

data2 = pd.DataFrame(val_2, columns=ns, index=ps)
ax = sns.heatmap(data2, cmap= 'coolwarm', vmin = val_2.min(), vmax = val_2.max())
plt.xlabel("Количество точек")
plt.ylabel("Степень кривой")
plt.savefig(f"plot_2_{p_start}_{p_end}_{n_start}_{n_end}")


