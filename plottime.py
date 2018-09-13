import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
#from collections import namedtuple


# Collect the data from the file, ignore empty lines
data = open('timeit.txt', 'r')

name = []
totalCPUtime = []

AllRuntime = dict()
FaceRuntime = dict()
CellRuntime = dict()

OrigRuntime = dict()
InPlaceRuntime = dict()
RegGridRuntime = dict()

for line in data:
    if 'CPU' in line:
        line = line.strip('\n')
        key,dummy,dummy,value,dummy = line.split()
        name.append(key)
        totalCPUtime.append(float(value))
        AllRuntime[key] = float(value)
        if 'face' in key:
           FaceRuntime[key] = float(value)
        if 'cell' in key:
           CellRuntime[key] = float(value)
        if 'orig' in key:
           OrigRuntime[key] = float(value)
        if 'inplace' in key:
           InPlaceRuntime[key] = float(value)
        if 'reggrid' in key:
           RegGridRuntime[key] = float(value)

fig, ax = plt.subplots()

bar_width = 0.35

y1 = OrigRuntime.values()
y2 = InPlaceRuntime.values()
y3 = RegGridRuntime.values()

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]
x3 = [x + bar_width for x in x2]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Orig')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='InPlace')
#plt.bar(x3, y3, width=bar_width, color='g', edgecolor='white', label='RegGrid')


ax.set_xlabel('Mesh Data Structure')
ax.set_ylabel('Runtime, secs')
ax.set_xticks(x1 + bar_width / 2)
ax.set_xticklabels(('Cells','Faces'))
ax.legend()

fig.tight_layout()
plt.show()
