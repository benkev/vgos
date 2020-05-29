'''
show_colormap.py
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#
# For indication of goodness of the multiple correlation coefficient,
# create rYlGn, a YlGn (Yellow-Green) colormap with reduced dynamic range
#
cmYlGn = plt.cm.get_cmap('YlGn')
rYlGn = ListedColormap(cmYlGn(np.linspace(0.0, 0.6, 256)), name='rYlGn')

#
# To plot the cross-correlation matrices, use inverted 'hot' colormap 
# with reduced dynamic range (From white throuhg yellow to dense red)
#
cmhot_r = plt.cm.get_cmap('hot_r')
hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')

x = np.arange(256)
y = np.zeros(256)
cols = plt.cm.jet(x)
iplot = 1

fig = plt.figure(figsize=(14,8))

for colmap in [rYlGn, hotr]:
    ax = plt.subplot(2, 1, iplot)

    for icol in range(256):
        ax.plot(x, y, marker='s', markersize=20, \
                # markerfacecolor=rYlGn.colors[icol]) 
                markerfacecolor=cols[icol]) 
    ax.set_xlim(0,255)
    ax.set_ylim(-1,+1)
    iplot += 1


fig.show()    
