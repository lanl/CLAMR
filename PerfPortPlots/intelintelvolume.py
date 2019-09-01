import matplotlib
import matplotlib.pyplot as plt
import numpy as np

xspace = [1, 3, 5, 7, 9, 11]
labels = ['cell', 'face', 'cell-\nin-place', 'face-\nin-place', 'regular-grid-\nby-cell', 'regular-grid-\nby-face']
FLOPs = [175, 40, 52, 36, 50, 29]
GBytes = [9, 46, 64, 48, 26, 45]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, FLOPs, width, label='FLOPs\n(in thousands)', color='cornflowerblue')
rects2 = ax.bar(x + width/2, GBytes, width, label='Memory volume\n(GBytes)', color='lightslategray')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Total Amount')
ax.set_title('Intel Skylake with Intel Compiler')
ax.set_xticks(x)
ax.set_xlim([-0.6, 6])
ax.set_xticklabels(labels, fontsize=14, rotation="vertical")
ax.legend(loc="upper right")
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
