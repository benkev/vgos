'''
verify_selbp.py
'''

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.signal import medfilt
# from matplotlib.colors import ListedColormap
# import phasecal
# from phasecal import mult_corr, write_xcorrmx, write_title, write_numbers
import os, sys, glob, copy
# import re
import datetime, time
# from functools import reduce
# import getopt

fields = [('year', 'i2'),  ('month', 'i2'), ('day', 'i2'), 
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('delay_s', 'f8'), ('source_name', 'S16'), ('scan_name', 'S16')]

dat_pcmt = np.loadtxt(
    '/data/geodesy/3713/pcc_datfiles/vt9360k2.pcmt.BCD.XY.dat', dtype=fields)
ndat = dat_pcmt.shape[0]

#
# Put time axis values in hours to t_hr
#
tyear =   dat_pcmt['year'].astype(int)
tmonth =  dat_pcmt['month'].astype(int)
tday =    dat_pcmt['day'].astype(int)
thour =   dat_pcmt['hour'].astype(int)
tminute = dat_pcmt['minute'].astype(int)
tsecond = dat_pcmt['second'].astype(int)

datim0 = datetime.datetime(tyear[0], tmonth[0], tday[0])
ttuple0 = datim0.timetuple()
tstamp0 = time.mktime(ttuple0)

# exp_doy_time = datim0.strftime('%j, %H:%M:%S')


t_sec = np.asarray(tsecond, dtype=float)

for itim in range(ndat):
    datim = datetime.datetime(tyear[itim], tmonth[itim], tday[itim], \
                              thour[itim], tminute[itim], tsecond[itim])
    ttuple = datim.timetuple()
    tstamp = time.mktime(ttuple)
    t_sec[itim] = tstamp

t_hr = (t_sec - tstamp0)/3600.    # Time in hours  




dels = dat_pcmt['delay_s']

plt.figure()
plt.plot(t_hr, dels, 'b.')
plt.grid(1)
plt.show()








