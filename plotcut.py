import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots()

x1, y1 = np.loadtxt('cell.cut2000', delimiter=',', unpack=True)
x2, y2 = np.loadtxt('face.cut2000', delimiter=',', unpack=True)
x3, y3 = np.loadtxt('cellinplace.cut2000', delimiter=',', unpack=True)
x4, y4 = np.loadtxt('faceinplace.cut2000', delimiter=',', unpack=True)
x5, y5 = np.loadtxt('reggridbycell.cut2000', delimiter=',', unpack=True)
x6, y6 = np.loadtxt('reggridbyface.cut2000', delimiter=',', unpack=True)

plt.plot(x1,y1, label='Original AMR Cell', color='black')
plt.plot(x2,y2, label='Original AMR Face', color='grey')
plt.plot(x3,y3, label='In-Place Cell', color='red')
plt.plot(x4,y4, label='In-Place Face', color='tomato')
plt.plot(x5,y5, label='Regular Grid Cell', color='blue')
plt.plot(x6,y6, label='Regular Grid Face', color='royalblue')

axes = plt.gca() # get current axes
axes.set_ylim([0,110])
#axes.set_xlim([0,70])

plt.xlabel('X',fontsize=18)
plt.ylabel('Height',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
plt.legend()
plt.setp(ax.get_legend().get_texts(), fontsize='14') # for legend text
fig.tight_layout()
plt.savefig("cutplanes.pdf")
plt.show()
