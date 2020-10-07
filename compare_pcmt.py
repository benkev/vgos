'''
compare_pcmt.py
'''
import numpy as np
from scipy.signal import medfilt
import matplotlib.pyplot as plt

fields = [('year', 'i4'),  ('month', 'i2'),  ('day', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('delay_s', 'f8'), ('source', 'a9'), ('delay_ps', 'a9')]

fname1 = '/data/geodesy/3618/pcc_datfiles/vt7254gs.pcmt.BCD.XY.dat'
fname2 = 'pcmt/vt7254gs.pcmt.BCD.XY.dat'

dat1 = np.loadtxt(fname1, dtype=fields)
dat2 = np.loadtxt(fname2, dtype=fields)

dl1 = dat1['delay_s']*1e12
dl2 = dat2['delay_s']*1e12


plt.figure(); plt.plot(dl1, 'b.'); plt.plot(dl2, 'g.'); plt.grid(1)







plt.show()







