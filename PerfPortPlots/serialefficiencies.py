import matplotlib
import matplotlib.pyplot as plt
import numpy as np

labels = ['cell', 'face', 'cell-\nin-place', 'face-\nin-place', 'regular-grid-\nby-cell', 'regular-grid-\nby-face']
Runtime = []
MBytes = [2, 19, 17, 22, 8, 20]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, Runtime, width, label='Runtime', color='cornflowerblue')
rects2 = ax.bar(x + width/2, MBytes, width, label='Bandwidth', color='lightslategray')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Efficiency (% of best)')
ax.set_title('Method Efficiencies for Serial Methods')
ax.set_xticks(x)
ax.set_xlim([-0.6, 6])
ax.set_xticklabels(labels, fontsize=14, rotation="vertical")
ax.legend(loc='upper right', prop={'size': 12})
fig.subplots_adjust(right=2.0)

"""
def autolabel(rects):
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
"""
fig.tight_layout()

plt.show()
