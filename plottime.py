#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
#from collections import namedtuple


# Collect the data from the file, ignore empty lines
data = open('timeit.txt', 'r')

name = []
totalCPUtime = []

Allruntime = dict()
Faceruntime = dict()
Cellruntime = dict()

for line in data:
    if 'CPU' in line:
        line = line.strip('\n')
        #print line
        key,dummy,dummy,value,dummy = line.split()
        name.append(key)
        totalCPUtime.append(value)
        #print name
        #print totalCPUtime
        Allruntime[key] = value
        if 'face' in key:
           Faceruntime[key] = value
        if 'cell' in key:
           Cellruntime[key] = value

print list(Allruntime)
print list(Faceruntime)
        #lines = [line.strip().split(': ') for line in f if len(line) > 1]



n_groups = 5

#means_men = (20, 35, 30, 35, 27)
#std_men = (2, 3, 4, 1, 2)
#
#means_women = (25, 32, 34, 20, 25)
#std_women = (3, 5, 2, 3, 3)

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
#error_config = {'ecolor': '0.3'}

rects1 = ax.bar(index, Allruntime.values(),  bar_width,
                alpha=opacity, color='b',
                label='Total CPU Time')
#
#rects2 = ax.bar(index + bar_width, means_women, bar_width,
#                alpha=opacity, color='r',
#                yerr=std_women, error_kw=error_config,
#                label='Women')
#
ax.set_xlabel('Mesh Data Structure')
ax.set_ylabel('Runtime, secs')
#ax.set_title('Scores by group and gender')
#ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(Allruntime.keys())
ax.legend()
#
#fig.tight_layout()
#plt.show()
