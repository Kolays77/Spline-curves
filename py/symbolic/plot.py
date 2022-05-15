import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rc('text', usetex = True)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

plt.rcParams["savefig.dpi"] = 300
plt.rcParams["savefig.bbox"] = 'tight'

x = [1, 2, 3, 4]
y = [3/6, 4/6, 5/6, 6/6]

plt.scatter(x, y, c="g")
plt.plot(x, y)

plt.yticks(y,
           [r"$B_{1p} w_1$", r"$B_{2p} w_2$", r"$B_{3p} w_3$", r"$B_{4p} w_4$"])
plt.show()
