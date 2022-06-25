import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
import pandas as pd
import sys
from matplotlib.colors import LogNorm, Normalize

# Использование: python3 plot_heat_map.py

plt.rcParams["savefig.bbox"] = 'tight'
plt.rcParams["savefig.dpi"] = 300

# Метод отрисовки температурной карты F(p, n), 
# где степень p меняется в пределах от p_start до p_end,
# число сегментов n  меняется в пределах от n_start до n_end.
# Параметры: 
# path - путь к файлу.
# файл структура:
# Первая строка: p_start p_end n_start n_end.
def plot_map(path):
    p_start = 0 
    p_end = 0
    n_start = 0
    n_end = 0

    values = []
    with open(path) as file:
        params = list(map(int, file.readline().split()))
        p_start = params[0]
        p_end = params[1]
        n_start = params[2]
        n_end = params[3]
        for line in file:
            l = list(map(float, line.split()))
            values.append(l)

    values = np.array(values)
    ps = list(range(p_start, p_end + 1))
    ns = list(range(n_start, n_end + 1))

    data1 = pd.DataFrame(values, columns=ns, index=ps)
    ax = sns.heatmap(data1, norm=LogNorm(), cmap= 'coolwarm', vmin = values.min(), vmax = values.max())
    plt.xlabel("Число сегментов сплайна", fontsize="x-large")
    plt.ylabel("Степень кривой", fontsize="x-large")
    plt.savefig(path)
    plt.close()

if __name__ == "__main__":
    path = sys.argv[1]
    plot_map(path)