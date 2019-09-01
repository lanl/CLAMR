"""import matplotlib
import matplotlib.pyplot as plt
import numpy as np

labels = ['cell', 'face', 'cell-\nin-place', 'face-\nin-place', 'regular-grid-\nby-cell', 'regular-grid-\nby-face']
OpenMP = [50, 95, 47, 100, 72, 83]
OpenCL = [9, 88, 62, 100, 38, 80]
MPI = [8]

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars
r1 = x
r2 = 

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/3, OpenMP, width, label='OpenMP', color='cornflowerblue')
rects2 = ax.bar(x, OpenCL, width, label='OpenCL', color='lightslategray')
rects2 = ax.bar(x + width/3, OpenCL, width, label='MPI', color='lightslategray')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Efficiency (% of best)')
ax.set_title('Method Efficiencies for Parallel Methods')
ax.set_xticks(x)
ax.set_xlim([-0.6, 6])
ax.set_xticklabels(labels, fontsize=14, rotation="vertical")
ax.legend(loc='upper right', prop={'size': 12})
fig.subplots_adjust(right=2.0)

"""
"""def autolabel(rects):
        Attach a text label above each bar in *rects*, displaying its height.
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

plt.show()
"""
# libraries
import numpy as np
import matplotlib.pyplot as plt
 
# set width of bar
barWidth = 0.25
labels = ['cell', 'face', 'cell\nin-place', 'face\nin-place', 'regular-grid\nby-face', 'regular-grid\nby-cell']
  
# set height of bar
bars1 = [1, 14, 13, 15, 14, 6]
bars2 = [1, 9, 12, 13, 0, 0]
bars3 = [3, 8, 13, 15, 0, 0]
   
# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
    
# Make the plot
plt.bar(r1, bars1, color='cornflowerblue', width=barWidth, edgecolor='white', label='Serial')
plt.bar(r2, bars2, color='forestgreen', width=barWidth, edgecolor='white', label='Heterogeneous\n(Serial+GPU)')
plt.bar(r3, bars3, color='slategrey', width=barWidth, edgecolor='white', label='Parallel\n(OMP+GPU)')
     
# Add xticks on the middle of the group bars
#plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(bars1))], labels)
plt.ylabel('Efficiency (% of best)')
plt.title('Bandwidth Architecture Efficiencies')
      
# Create legend & Show graphic
plt.legend(prop={'size' : 10})
plt.show()
