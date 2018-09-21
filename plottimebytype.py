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

AllMeshtime = dict()
FaceMeshtime = dict()
CellMeshtime = dict()

OrigMeshtime = dict()
InPlaceMeshtime = dict()
RegGridMeshtime = dict()

for line in data:
    #if 'CPU' in line:
    if 'finite_difference' in line:
        line = line.strip('\n')
        #key,dummy,dummy,value,dummy = line.split()
        key,dummy,value,dummy = line.split()
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
    if 'mesh_timer_bidir  ' in line:
        line = line.strip('\n')
        #key,dummy,dummy,value,dummy = line.split()
        key,dummy,value,dummy = line.split()
        name.append(key)
        if 'face' in key:
           FaceMeshtime[key] = float(value)
        if 'cell' in key:
           CellMeshtime[key] = float(value)
        if 'orig' in key:
           OrigMeshtime[key] = float(value)
        if 'inplace' in key:
           InPlaceMeshtime[key] = float(value)
        if 'reggrid' in key:
           RegGridMeshtime[key] = float(value)

fig, ax = plt.subplots()

bar_width = 0.35

y1 = [ CellRuntime["origcell"]/CellRuntime["origcell"],CellRuntime["cellinplace"]/CellRuntime["origcell"],CellRuntime["reggridbycell"]/CellRuntime["origcell"] ]
y2 = [ FaceRuntime["origface"]/CellRuntime["origcell"],FaceRuntime["faceinplace"]/CellRuntime["origcell"],FaceRuntime["reggridbyfaces"]/CellRuntime["origcell"] ]
y3 = [ CellMeshtime["origcell"]/CellRuntime["origcell"],CellMeshtime["cellinplace"]/CellRuntime["origcell"],CellMeshtime["reggridbycell"]/CellRuntime["origcell"] ]
y4 = [ FaceMeshtime["origface"]/CellRuntime["origcell"],FaceMeshtime["faceinplace"]/CellRuntime["origcell"],FaceMeshtime["reggridbyfaces"]/CellRuntime["origcell"] ]

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Cells')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='Faces')
plt.bar(x1, y3, width=bar_width, color='k', edgecolor='white', label='Mesh')
plt.bar(x2, y4, width=bar_width, color='k', edgecolor='white')
#plt.bar(x3, y3, width=bar_width, color='g', edgecolor='white', label='RegGrid')

axes = plt.gca() # get current axes
axes.set_ylim([0,1.5])

ax.set_xlabel('Mesh Data Structure',fontsize=18)
ax.set_ylabel('Runtime relative to Original Cell AMR',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
ax.set_xticks(x1 + bar_width)
ax.set_xticklabels(('Original AMR','InPlace','RegGrid'))
ax.legend(loc=2)

plt.setp(ax.get_legend().get_texts(), fontsize='16') # for legend text

fig.tight_layout()
plt.savefig("plottimebytype.pdf")
plt.show()
