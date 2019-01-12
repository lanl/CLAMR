import numpy as np
import matplotlib.pyplot as plt

# Collect the data from the file, ignore empty lines
data = open('linecount.txt', 'r')

name = []
totalCPUtime = []

AllRuntime = dict()
FaceLineCount = dict()
CellLineCount = dict()

OrigLineCount = dict()
InPlaceLineCount = dict()
RegGridLineCount = dict()

for line in data:
    if 'State' in line:
       continue
    if '^calc' in line:
       continue
    if 'finite_difference' in line:
        line = line.strip('\n')
        #print(len(line.split()))
        if len(line.split()) == 4:
           string,dummy,dummy,value = line.split()
           if 'face' in string:
              if 'via_faces' in string:
                 FaceLineCount["origface"] = float(value)
              if 'face_in_place' in string:
                 FaceLineCount["faceinplace"] = float(value)
              if 'regular_cells_by_faces' in string:
                 FaceLineCount["reggridbyfaces"] = float(value)
           else:
              if 'difference.cpp' in string:
                 CellLineCount["origcell"] = float(value)
              if 'cell_in_place' in string:
                 CellLineCount["cellinplace"] = float(value)
              if 'regular_cells.cpp' in string:
                 CellLineCount["reggridbycell"] = float(value)


fig, ax = plt.subplots()

bar_width = 0.35

#y1 = [ CellLineCount["origcell"]/CellLineCount["origcell"],CellLineCount["cellinplace"]/CellLineCount["origcell"],0 ]
#y2 = [ FaceLineCount["origface"]/CellLineCount["origcell"],FaceLineCount["faceinplace"]/CellLineCount["origcell"],0 ]

#y1 = [ CellLineCount["origcell"]/CellLineCount["origcell"],CellLineCount["cellinplace"]/CellLineCount["origcell"],CellLineCount["reggridbycell"]/CellLineCount["origcell"] ]
#y2 = [ FaceLineCount["origface"]/CellLineCount["origcell"],FaceLineCount["faceinplace"]/CellLineCount["origcell"],FaceLineCount["reggridbyfaces"]/CellLineCount["origcell"] ]
y1 = [ CellLineCount["origcell"],CellLineCount["cellinplace"],CellLineCount["reggridbycell"] ]
y2 = [ FaceLineCount["origface"],FaceLineCount["faceinplace"],FaceLineCount["reggridbyfaces"] ]

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Cells')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='Faces')

axes = plt.gca() # get current axes
axes.set_ylim([0,500])

ax.set_xlabel('Mesh Data Structure',fontsize=18)
ax.set_ylabel('Linecount',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
ax.set_xticks(x1 + bar_width)
ax.set_xticklabels(('Original AMR','InPlace','RegGrid'))
ax.legend(loc=1)

plt.setp(ax.get_legend().get_texts(), fontsize='16') # for legend text

fig.tight_layout()
plt.savefig("plotLineCount.pdf")
plt.savefig("plotLineCount.svg")
plt.show()
