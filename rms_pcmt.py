'''
rms_pcmt.py

Comparison of the averaged cable delays read from the PCMT Files.

The comparison is based on the RMS of their difference.

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import sys
import datetime, time

fields = [('year', 'i4'),  ('month', 'i2'),  ('day', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('delay_s', 'f8'), ('source', 'a9'), ('delay_ps', 'a9')]


pcmt1 = sys.argv[1]   # First pcmt file name
pcmt2 = sys.argv[2]   # Second pcmt file name

dat1 = np.loadtxt(pcmt1, dtype=fields)
dat2 = np.loadtxt(pcmt2, dtype=fields)

ndat = len(dat1)
tyear =   dat1['year'].astype(int)
tmonth =  dat1['month'].astype(int)
tday =    dat1['day'].astype(int)
thour =   dat1['hour'].astype(int)
tminute = dat1['minute'].astype(int)
tsecond = dat1['second'].astype(int)

datim0 = datetime.datetime(tyear[0], tmonth[0], tday[0], thour[0], \
                           tminute[0], tsecond[0])
ttuple0 = datim0.timetuple()
tstamp0 = time.mktime(ttuple0)

exp_doy_time = datim0.strftime('%j, %H:%M:%S')

t_sec = np.asarray(tsecond, dtype=float)

for itim in range(ndat):
    datim = datetime.datetime(tyear[itim], tmonth[itim], tday[itim], \
                               thour[itim], tminute[itim], tsecond[itim])
    ttuple = datim.timetuple()
    tstamp = time.mktime(ttuple)
    t_sec[itim] = tstamp
    
t_hr = (t_sec - tstamp0)/3600.    # Time in hours  




dl1 = dat1['delay_s']*1e12
dl2 = dat2['delay_s']*1e12

#
# Remove trend using the median filter
#
delt_tr = dl1 - dl2             # Difference with trend
delt_mf = medfilt(delt_tr, 21)  # Median filtered
delt = delt_tr - delt_mf        # Trend (median) removed

rms = np.sqrt(np.mean(np.square(delt)))

rms_mf = medfilt(rms, 21)

print('rms: ' + '%8.3f\n' % rms)


fig = plt.figure(figsize=(12,10))
#fig = plt.figure(figsize=(16,5))
fig.suptitle('Averaged Cable Delays Comparison (PCMT Files). RMS = ' \
             + '%6.3f\n' % rms, fontsize=16)

ax1 = plt.subplot(221)
ax1.plot(t_hr, dl1, 'b.'); ax1.grid(1)
plt.ylabel('Cable 1 Delay (ps)')
plt.xlabel('Time (hours)')

ax2 = plt.subplot(222)
ax2.plot(t_hr, dl2, 'g.'); ax2.grid(1)
plt.ylabel('Cable 2 Delay (ps)')
plt.xlabel('Time (hours)')

ax3 = plt.subplot(223)
ax3.plot(t_hr, delt_tr, 'k.'); ax3.grid(1)
ax3.plot(t_hr, delt_mf, 'r', lw=3)
plt.ylabel('Difference (ps)')
plt.xlabel('Time (hours)')

ax4 = plt.subplot(224)
ax4.plot(t_hr, delt, 'r.'); ax4.grid(1)
ax4.set_ylim(2*delt.min(), 2*delt.max())  # Vertically narrow the error plot 
plt.ylabel('Difference untrended (ps)')
plt.xlabel('Time (hours)')

fig.show()





