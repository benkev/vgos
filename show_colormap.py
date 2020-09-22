'''
show_colormap.py
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

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
cmap = ListedColormap([
    '#08306b', '#08519c', '#2171b5', '#4292c6',
    '#6baed6', '#9ecae1', '#c6dbef', '#deebf7',
    '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a',
    '#ef3b2c', '#cb181d', '#a50f15', '#67000d',
])

cmap1 = ListedColormap([
    '#08306b', '#08519c', '#2171b5', '#4292c6',
    '#6baed6', '#9ecae1', '#c6dbef', '#deebf7',
    '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a',
    '#ef3b2c', '#cb181d', '#a50f15'
])

cmjet = plt.cm.get_cmap('jet')
lcmjet = ListedColormap(cmjet(np.linspace(0.1, 0.9, 16)), name='lcmjet')

cmhsv = plt.cm.get_cmap('hsv')
lcmhsv = ListedColormap(cmhsv(np.linspace(0., 0.35, 8)), name='lcmhsv')


r = np.linspace(1, 0, 16)
g = np.linspace(0, 1, 16)
gr = np.zeros((16,4)) 
gr[:,0] = r
gr[:,1] = g
gr[:,2] = 0 
gr[:,3] = 1

lcmgr = ListedColormap(gr)


iplot = 1

# fig = plt.figure(figsize=(14,8))
fig = plt.figure()


for colmap in [lcmjet, hotr, lcmhsv, lcmgr]:
    ax = fig.add_subplot(4, 1, iplot, xlim=(0,255), ylim=(0,1))
    #ax.set_axis_off()
    col = PatchCollection([Rectangle((x, 0), 1, 1) for x in range(256)])
    col.set_array(np.arange(256))
    col.set_cmap(colmap)
    ax.add_collection(col)

    iplot += 1


fig.show()    


