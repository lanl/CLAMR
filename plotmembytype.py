import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
#from collections import namedtuple


# Collect the data from the file, ignore empty lines
data = open('timeit.txt', 'r')

name = []
totalCPUtime = []

AllMemused = dict()
FaceMemused = dict()
CellMemused = dict()

OrigMemused = dict()
InPlaceMemused = dict()
RegGridMemused = dict()

for line in data:
    #if 'CPU' in line:
    if 'Memory_used' in line:
        line = line.strip('\n')
        #key,dummy,dummy,value,dummy = line.split()
        key,dummy,value,dummy = line.split()
        name.append(key)
        totalCPUtime.append(float(value))
        AllMemused[key] = float(value)
        if 'face' in key:
           FaceMemused[key] = float(value)
        if 'cell' in key:
           CellMemused[key] = float(value)
        if 'orig' in key:
           OrigMemused[key] = float(value)
        if 'inplace' in key:
           InPlaceMemused[key] = float(value)
        if 'reggrid' in key:
           RegGridMemused[key] = float(value)

fig, ax = plt.subplots()

bar_width = 0.35

y1 = [ CellMemused["origcell"]/CellMemused["origcell"],CellMemused["cellinplace"]/CellMemused["origcell"],CellMemused["reggridbycell"]/CellMemused["origcell"] ]
y2 = [ FaceMemused["origface"]/CellMemused["origcell"],FaceMemused["faceinplace"]/CellMemused["origcell"],0 ]
#y2 = [ FaceMemused["origface"]/CellMemused["origcell"],FaceMemused["faceinplace"]/CellMemused["origcell"],FaceMemused["reggridbyfaces"]/CellMemused["origcell"] ]
#y1 = [ CellMemused["origcell"]/CellMemused["origcell"],CellMemused["cellinplace"]/CellMemused["origcell"] ]
#y2 = [ FaceMemused["origface"]/CellMemused["origcell"],FaceMemused["faceinplace"]/CellMemused["origcell"] ]

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Cells')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='Faces')
#plt.bar(x3, y3, width=bar_width, color='g', edgecolor='white', label='RegGrid')

axes = plt.gca() # get current axes
#axes.set_ylim([0,1.2])

ax.set_xlabel('Mesh Data Structure',fontsize=18)
ax.set_ylabel('Memory used relative to Original Cell AMR',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
ax.set_xticks(x1 + bar_width)
ax.set_xticklabels(('Original AMR','InPlace','RegGrid'))
ax.legend(loc=1)

plt.setp(ax.get_legend().get_texts(), fontsize='16') # for legend text

fig.tight_layout()
plt.savefig("plotmembytype.pdf")
plt.show()
