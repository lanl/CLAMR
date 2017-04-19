#!/usr/bin/env python

import os
import sys
import string
import numpy
import matplotlib.pyplot as plt
#axes.titlesize : large
#axes.labelsize : medium
#font = {'family' : 'normal',
#    'weight' : 'bold',
#        'size'   : 67}
# Must define the files:
if len(sys.argv) == 1:
    print "Syntax: python2.7 self_plot.py <FILES>"; exit(1)

# Initialize Plot:
fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)

# Loop through the files:
for tblfile in sys.argv[1:]:

    # Initialize empty lists:
    #nodes   = []
    #ranks   = []
    threads = []
    speed_omp   = []
    speed_homp  = []
    ideal_speed     = []
    
    # Load file into an array:
    all_lines = numpy.loadtxt(tblfile, comments="#")

    # Loop through each line (list of items) in the "all_lines" array:
    for line in all_lines:

        # Append new info to lists
        #nodes.append(int(line[0]))
        #ranks.append(int(line[1]))
        threads.append(int(line[0]))
        speed_omp.append(float(line[1]))
        speed_homp.append(float(line[2]))
        ideal_speed.append(int(line[3]))

    #for i in range(1,len(speed_omp)):
    for i in range(1,5):
       speed_omp[i] = speed_omp[0]/speed_omp[i]
       speed_homp[i] = speed_homp[0]/speed_homp[i]

    speed_omp[0] = 1.0
    speed_homp[0] = 1.0

    #err.append(0.0)
    
    ax = axs #ax = axs[0,0]

    #ax.errorbar(nodes, memuse, err, fmt='s--', label=curve_label)
   # plt.text(2.3,3, 'Up to 89% parallel efficiency', fontsize=35)
    plt.plot(threads, ideal_speed, 's--', linewidth= 3.5, label = 'Ideal Speed-Up')
    plt.plot(threads, speed_homp, 's--', linewidth= 3.5,label = 'High-level OpenMP')
    plt.plot(threads, speed_omp, 's--', linewidth= 3.5,label = 'MPI')



# Shrink current axis by 20%
if False:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
else:
    ax.legend(loc='best',fontsize=45)

ax.set_title('High-level OpenMP VS MPI on CLAMR',fontsize=50)
ax.set_xlabel('Number of Threads or MPI Processes(Respectively)',fontsize=50)
ax.set_ylabel('Speed Up',fontsize=50)

for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(40)
for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(40)

#ax.set_xlim([400,1000])
#ax.set_ylim([10**-1,10**6])
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.show()

#nthreads omp_sup hlom_sup ideal_sup
