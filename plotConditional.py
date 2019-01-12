import numpy as np
import matplotlib.pyplot as plt

# Collect the data from the file, ignore empty lines
data = open('linecount.txt', 'r')

name = []
totalCPUtime = []

AllRuntime = dict()
FaceConditional = dict()
CellConditional = dict()

OrigConditional = dict()
InPlaceConditional = dict()
RegGridConditional = dict()

for line in data:
    if 'State' in line:
       continue
    if '^calc' in line:
       continue
    if 'finite_difference' in line:
        line = line.strip('\n')
        #print(len(line.split()))
        if len(line.split()) == 2:
           string,value = line.split()
           if 'face' in string:
              if 'via_faces' in string:
                 FaceConditional["origface"] = float(value)
              if 'face_in_place' in string:
                 FaceConditional["faceinplace"] = float(value)
              if 'regular_cells_by_faces' in string:
                 FaceConditional["reggridbyfaces"] = float(value)
           else:
              if 'difference.cpp' in string:
                 CellConditional["origcell"] = float(value)
              if 'cell_in_place' in string:
                 CellConditional["cellinplace"] = float(value)
              if 'regular_cells.cpp' in string:
                 CellConditional["reggridbycell"] = float(value)


fig, ax = plt.subplots()

bar_width = 0.35

#y1 = [ CellConditional["origcell"]/CellConditional["origcell"],CellConditional["cellinplace"]/CellConditional["origcell"],0 ]
#y2 = [ FaceConditional["origface"]/CellConditional["origcell"],FaceConditional["faceinplace"]/CellConditional["origcell"],0 ]

#y1 = [ CellConditional["origcell"]/CellConditional["origcell"],CellConditional["cellinplace"]/CellConditional["origcell"],CellConditional["reggridbycell"]/CellConditional["origcell"] ]
#y2 = [ FaceConditional["origface"]/CellConditional["origcell"],FaceConditional["faceinplace"]/CellConditional["origcell"],FaceConditional["reggridbyfaces"]/CellConditional["origcell"] ]
y1 = [ CellConditional["origcell"],CellConditional["cellinplace"],CellConditional["reggridbycell"] ]
y2 = [ FaceConditional["origface"],FaceConditional["faceinplace"],FaceConditional["reggridbyfaces"] ]

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Cells')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='Faces')

axes = plt.gca() # get current axes
#axes.set_ylim([0,500])

ax.set_xlabel('Mesh Data Structure',fontsize=18)
ax.set_ylabel('Conditionals',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
ax.set_xticks(x1 + bar_width)
ax.set_xticklabels(('Original AMR','InPlace','RegGrid'))
ax.legend(loc=1)

plt.setp(ax.get_legend().get_texts(), fontsize='16') # for legend text

fig.tight_layout()
plt.savefig("plotConditional.pdf")
plt.savefig("plotConditional.svg")
plt.show()
