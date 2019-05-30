import os
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from pylab import rcParams
from sklearn.cluster import dbscan
from tefingerprint.cluster import DBICAN, SDBICAN

# directory of this script
directory = os.path.dirname(os.path.abspath(__file__))

# figure 1: comparison of DBSCAN DBSCAN* and DBICAN

# read tip positions
points = np.array([1, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 8, 11, 12,
                   12, 13, 13, 13, 14, 14, 14, 15, 15, 16, 18,
                   23, 25, 26, 27, 27, 28, 29, 30, 31, 33])

# variables
mpts = 10
eps = 5

# DBSCAN labels
cores, dbscan_labels = dbscan(points.reshape(-1, 1),
                              eps=eps,
                              min_samples=mpts)

# DBSCAN* labels are same for core points but otherwise -1
dbscanx_labels = np.zeros(len(points), dtype=np.int) - 1
dbscanx_labels[cores] = dbscan_labels[cores]

# DBICAN labels
dbican = DBICAN(min_points=mpts, epsilon=eps)
dbican.fit(points)
dbican_labels = dbican.labels()

# plotting
labels_list = [dbscan_labels, dbscanx_labels, dbican_labels]
name_list = [r'DBSCAN', r'DBSCAN*', 'DBICAN']
title = r'Comparison of Algorithms with $m_{pts}=10$ and $\varepsilon=5$'
legend_labels = ['.', '1', '2']
x_max_ofset = 2
height = 2.5
width = 6
n_row = 3
n_col = 1

rcParams['figure.figsize'] = width, height
x_max = np.max(points) + x_max_ofset
position = np.arange(x_max)
colour = np.zeros(position.shape, dtype='U2')
grey = '#C8C8C8'
colour[:] = grey

fig, ax = plt.subplots(n_row, n_col, sharex='col', sharey='row')

for row in range(n_row):
    name = name_list[row]
    labels = labels_list[row]

    for l in [-1, 0, 1]:
        if l is -1:
            colour = grey
        else:
            colour = 'C{0}'.format(l)
        counts = np.zeros(x_max, dtype=np.int)
        counts_ = Counter(points[labels == l])
        for i in range(len(position)):
            if i in counts_:
                counts[i] = counts_[i]

        ax[row].bar(position,
                    counts,
                    width=1,
                    color=colour,
                    tick_label=position)
        if row is 0:
            ax[row].arrow(23, 2.5, dx=0, dy=-1,
                          length_includes_head=True,
                          head_width=0.5,
                          head_length=0.2)

        if row is 0:
            ax[row].set_title(title)
        ax[row].get_xaxis().set_ticks([])
        ax[row].get_yaxis().set_ticks([])
        ax[row].set_ylabel(name)

plt.savefig(directory+'/DBICAN.png')


# figure 2: comparison of DBICAN and SDBICAN

# read tip positions
points = np.array([1, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 8, 11, 12,
                   12, 13, 13, 13, 14, 14, 14, 15, 15, 16, 18])

# variables
mpts = 10

# plot
legend_labels = ['.', '1', '2']
x_max_ofset = 2
height = 5
width = 6
n_row = 10
n_col = 2

rcParams['figure.figsize'] = width, height
x_max = np.max(points) + x_max_ofset
position = np.arange(x_max)
colour = np.zeros(position.shape, dtype='U2')
grey = '#C8C8C8'
colour[:] = grey

fig, ax = plt.subplots(n_row, n_col, sharex='col', sharey='row')

for col in range(n_col):

    for row in range(n_row):
        mpts = 10
        eps = row + 1
        if col is 0:
            model = DBICAN(mpts, eps)
        elif col is 1:
            model = SDBICAN(mpts, eps)

        model.fit(points)
        labels = model.labels()

        for l in [-1, 0, 1]:
            if l is -1:
                colour = grey
            else:
                colour = 'C{0}'.format(l)
            counts = np.zeros(x_max, dtype=np.int)
            counts_ = Counter(points[labels == l])
            for i in range(len(position)):
                if i in counts_:
                    counts[i] = counts_[i]

            ax[row, col].bar(position,
                             counts,
                             width=1,
                             color=colour,
                             tick_label=position)
            ax[row, col].get_xaxis().set_ticks([])
            ax[row, col].get_yaxis().set_ticks([])

            if row is 0:
                ax[row, col].set_title('SDBICAN' if col else 'DBICAN')
            if col is 0:
                ax[row, col].set_ylabel(r'$\varepsilon={0}$'.format(eps))

plt.savefig(directory+'/SDBICAN.png')
