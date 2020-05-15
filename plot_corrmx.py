'''
plot_corrmx.py

Plot matrix of Pearson's correlation coefficients for all channels.
Usage:

plot_corrmx.py <4-digit experiment code> <experiment name, VT and 4 digits> \
                    <one-letter station name>
Like:
plot_corrmx.py 3638 VT8078 Y
'''


import numpy as np
from scipy.signal import medfilt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os, sys, glob, copy

exprcode = str(sys.argv[1])           # Like 3694
exprname = str(sys.argv[2]).upper()   # Like VT9162
station =  sys.argv[3].upper()        # A letter, either of E, G, H, I, V, T 


fields = [('year', 'i2'),  ('doy', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('phase_midband', 'f8'), ('phase_dc', 'f8'),
          ('delay_ps', 'f8'), ('phase_rmse', 'f8'),
          ('scan', 'S16'), ('source', 'S16'), ('station', 'S2'),
          ('azimuth', 'f8'), ('elevation', 'f8')]


ddir = '/data/geodesy/' + exprcode + '/pcc_datfiles'
if os.path.isdir(ddir):
    datdir = copy.copy(ddir)
elif os.path.isdir(ddir + '_jb'):
    datdir = ddir + '_jb'
else:
    print('Neither "'+datdir+'" nor "'+datdir+ '_jb" exists. Skipped.')
    sys.exit(1) # ============================================== >>>

fnames = glob.glob(datdir + '/bandmodel.??????.' + station + '.?.?.dat')
fnames.sort()

nfile = len(fnames)

datlist = []
for ix in range(nfile):
    datlist.append(np.loadtxt(fnames[ix], dtype=fields))

dat = np.array(datlist)  # Convert list of struct. arrays to array
ndat = np.size(dat,1)
delps = np.zeros((nfile,ndat), dtype=float)  # Cable delays (ps)
for ix in range(nfile):
    delps[ix,:] = dat[ix]['delay_ps']

nbandpol = nfile # Assume there are 8 files for each band and pol


#
# Compute 8x8 correlation matrix of rows of dps[8,878]
#
Rxx_full = np.corrcoef(delps)
for ix in range(nbandpol):
	Rxx_full[ix,:(ix+1)] = np.NaN

# Rxx = np.delete(np.delete(Rxx_full, iy, axis=0), iy, axis=1)

#
# Use inverted 'hot' colormap with frduced dynamic range
# (From white throuhg yellow to dense red)
#
cmhot_r = plt.cm.get_cmap('hot_r')
hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')


plt.figure();
plt.imshow(Rxx_full, interpolation='none', cmap=hotr);
#plt.imshow(Rxx_full, interpolation='none', cmap=plt.cm.rainbow);
#plt.pcolormesh(Rxx_full, cmap=plt.cm.jet, offset_position='data');
plt.xticks(np.arange(8), ['AX', 'AY', 'BX', 'BY', 'CX', 'CY', 'DX', 'DY'])
plt.tick_params(axis='x', labeltop='on')
plt.yticks(np.arange(8), ['AX', 'AY', 'BX', 'BY', 'CX', 'CY', 'DX', 'DY'])
plt.grid(1)
plt.colorbar(shrink=0.8)
# plt.title('Pearson''s Correlations. Code ' + exprcode + 
#		  ' Experiment ' + exprname + '  Station ' + station)
plt.figtext(0.1, 0.95, 'Pearson\'s Correlations. Station ' + station +
            ', Code ' + exprcode + ', Experiment ' + exprname)
plt.show()

sys.exit(0)



plt.figure()
plt.plot(y, 'ko', ms=3, alpha=0.7, label='X Original'); plt.grid(1)
plt.plot(yf, 'r', lw=2, alpha=0.7, label='X Median Filtered')
#plt.plot(py, 'gx', alpha=0.7, label='dX/dt')
plt.plot(pyf, 'b', ms=4, alpha=0.5, label='dXmf/dt')  
plt.title('Cable delay, ' + dname + '.Y')
plt.ylabel('delay (ps)'); plt.xlabel('time')
plt.legend()

plt.show()


xmean =  np.mean(x)
ymean =  np.mean(y)
x0 = x - xmean
y0 = y - ymean
covxy = np.mean(x0*y0)
sig_x = np.std(x)
sig_y = np.std(y)

cor1 = covxy/(sig_x*sig_y) # Pearson's correlation coefficient

cor2 = np.corrcoef(x, y)[1,0]

