import numpy as np
import matplotlib.pyplot as plt

# Collect the data from the file, ignore empty lines
data = open('linecount.txt', 'r')

name = []
totalCPUtime = []

AllRuntime = dict()
FaceCCN = dict()
CellCCN = dict()

OrigCCN = dict()
InPlaceCCN = dict()
RegGridCCN = dict()

for line in data:
    if 'State' in line:
       continue
    if '^calc' in line:
       continue
    if 'finite_difference' in line:
        line = line.strip('\n')
        #print(len(line.split()))
        if len(line.split()) == 6:
           dummy,dummy,value,dummy,dummy,string = line.split()
           if 'face' in string:
              if 'via_faces' in string:
                 FaceCCN["origface"] = float(value)
              if 'face_in_place' in string:
                 FaceCCN["faceinplace"] = float(value)
              if 'regular_cells_by_faces' in string:
                 FaceCCN["reggridbyfaces"] = float(value)
           else:
              print(line)
              print(string)
              if 'difference.cpp' in string:
                 CellCCN["origcell"] = float(value)
              if 'cell_in_place' in string:
                 CellCCN["cellinplace"] = float(value)
              if 'regular_cells.cpp' in string:
                 CellCCN["reggridbycell"] = float(value)


fig, ax = plt.subplots()

bar_width = 0.35

#y1 = [ CellCCN["origcell"]/CellCCN["origcell"],CellCCN["cellinplace"]/CellCCN["origcell"],0 ]
#y2 = [ FaceCCN["origface"]/CellCCN["origcell"],FaceCCN["faceinplace"]/CellCCN["origcell"],0 ]

#y1 = [ CellCCN["origcell"]/CellCCN["origcell"],CellCCN["cellinplace"]/CellCCN["origcell"],CellCCN["reggridbycell"]/CellCCN["origcell"] ]
#y2 = [ FaceCCN["origface"]/CellCCN["origcell"],FaceCCN["faceinplace"]/CellCCN["origcell"],FaceCCN["reggridbyfaces"]/CellCCN["origcell"] ]

print(CellCCN)

y1 = [ CellCCN["origcell"],CellCCN["cellinplace"],CellCCN["reggridbycell"] ]
y2 = [ FaceCCN["origface"],FaceCCN["faceinplace"],FaceCCN["reggridbyfaces"] ]

x1 = np.arange(len(y1))
x2 = [x + bar_width for x in x1]

plt.bar(x1, y1, width=bar_width, color='b', edgecolor='white', label='Cells')
plt.bar(x2, y2, width=bar_width, color='r', edgecolor='white', label='Faces')

axes = plt.gca() # get current axes
#axes.set_ylim([0,500])

ax.set_xlabel('Mesh Data Structure',fontsize=18)
ax.set_ylabel('Cyclomatic Complexity',fontsize=18)
ax.tick_params(axis = 'both', labelsize = 14)
ax.set_xticks(x1 + bar_width)
ax.set_xticklabels(('Original AMR','InPlace','RegGrid'))
ax.legend(loc=1)

plt.setp(ax.get_legend().get_texts(), fontsize='16') # for legend text

fig.tight_layout()
plt.savefig("plotCCN.pdf")
plt.show()
