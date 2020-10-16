'''
compare_pcmt.py
'''
import numpy as np
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import re, os, sys, glob

fields = [('year', 'i4'),  ('month', 'i2'),  ('day', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('delay_s', 'f8'), ('source', 'a9'), ('delay_ps', 'a9')]
st1to2 = {
    'E' : 'wf', # Westford
    'F' : 'w2', # Westford2 (defunct)
    'G' : 'gs', # GGAO12M (Goddard)
    'H' : 'k2', # Kokee
    'V' : 'ws', # Wettzell
    'Y' : 'yj', # Yebes
    'I' : 'is', # Ishioka
    'S' : 'oe', # Onsala-Northeast
    'T' : 'ow', # Onsala-Southwest
    'M' : 'mg'  # MGO (MacDonald)
}

st2to1 = {
    'wf' : 'E', # Westford
    'w2' : 'F', # Westford2 (defunct)
    'gs' : 'G', # GGAO12M (Goddard)
    'k2' : 'H', # Kokee
    'ws' : 'V', # Wettzell
    'yj' : 'Y', # Yebes
    'is' : 'I', # Ishioka
    'oe' : 'S', # Onsala-Northeast
    'ow' : 'T', # Onsala-Southwest
    'mg' : 'M'  # MGO (MacDonald)
}


pdname = sys.argv[1]   # Input directory name (containing pcmt files)

with open('pcmt.txt') as f: pcmts = f.readlines()

fout = open(pdname + 'pcmt_stat.txt', 'w')

patt =  r'(?:\/[A-y]{2}[0-9]{4}|\/[A-y][0-9]{5})'
patt += r'[a-y][0-9a-y]\.pcmt\.'     # [ABCD]{1,4}\.[XY]{1,2}\.dat'

patt1 = r'(?:\/[A-y]{2}[0-9]{4}|\/[A-y][0-9]{5})[a-y][0-9a-y]'

rmsl = []

for pcmt in pcmts:
    
    pcmt = pcmt[:-1]                  # Remove '\n' at the end
    pcmt2 = re.findall(patt, pcmt)[0]
    pcmt2 = pcmt2[3:]
    pcmt2 = glob.glob(pdname + '*' + pcmt2 + '*.dat')

    if len(pcmt2) == 0: # Ignore, if no such file in the input directory
        continue   # ====================================================== >>>
    
    pcmt2 = pcmt2[0]

    dat1 = np.loadtxt(pcmt,  dtype=fields)
    dat2 = np.loadtxt(pcmt2, dtype=fields)

    dl1 = dat1['delay_s']*1e12
    dl2 = dat2['delay_s']*1e12
    Nsample = len(dl1)
    
    rms = np.sqrt(np.mean(np.square(dl1 - dl2)))
    

    expst = re.findall(patt1, pcmt) # Experiment name and 2-char station code
    expst = expst[0][1:]
    fout.write(expst + '  ' + '%8.3f\n' % rms)
    rmsl.append(rms)

    print(pcmt2)

    #plt.figure(); plt.plot(dl1, 'b.'); plt.plot(dl2, 'g.'); plt.grid(1)

avg = np.mean(rmsl)
print('mean rms: ' + '%8.3f\n' % avg)
fout.write('mean rms: ' + '%8.3f\n' % avg)




#plt.show()


fout.close()




