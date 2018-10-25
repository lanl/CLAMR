import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots()

x1, y1 = np.loadtxt('cell.cut1000', delimiter=',', unpack=True)
x2, y2 = np.loadtxt('face.cut1000', delimiter=',', unpack=True)
x3, y3 = np.loadtxt('cellinplace.cut1000', delimiter=',', unpack=True)
x4, y4 = np.loadtxt('faceinplace.cut1000', delimiter=',', unpack=True)
x5, y5 = np.loadtxt('reggridbycell.cut1000', delimiter=',', unpack=True)
x6, y6 = np.loadtxt('reggridbyface.cut1000', delimiter=',', unpack=True)

plt.scatter(x1,y1, marker="+", label='Original AMR Cell', color='black')
plt.scatter(x2,y2, marker="x", label='Original AMR Face', color='grey')
plt.scatter(x3,y3, marker="1", label='In-Place Cell', color='red')
plt.scatter(x4,y4, marker="2", label='In-Place Face', color='tomato')
plt.scatter(x5,y5, marker="3", label='Regular Grid Cell', color='blue')
plt.scatter(x6,y6, marker="4", label='Regular Grid Face', color='royalblue')

axes = plt.gca() # get current axes
axes.set_ylim([0,90])
axes.set_xlim([0,15])

plt.xlabel('Radius',fontsize=18)
plt.ylabel('Height',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
plt.legend()
plt.setp(ax.get_legend().get_texts(), fontsize='14') # for legend text
fig.tight_layout()
plt.savefig("cutplanesall.pdf")
plt.savefig("cutplanesall.jpg")
plt.show()
