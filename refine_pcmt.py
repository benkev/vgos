'''
refine_pcmt.py

The pcmt1.txt file is created by the command:
$ ls /data/geodesy/????/pcc_datfiles*/*.pcmt.*.dat > pcmt1.dat

It may contain unneeded files. This program cleand the pcmt1.txt file 
of the garbage and saves the result in the pcmt.txt file.

'''
import numpy as np
import re


pcmt_raw = np.loadtxt('pcmt1.txt', dtype='a128', delimiter='\n')

patt =  r'\/data\/geodesy'
patt += r'\/[0-9]{4}\/pcc_(?:datfiles|datfiles_jb)\/'
patt += r'(?:[A-y]{2}[0-9]{4}|[A-y][0-9]{5})'
patt += r'[0-9A-y]{2}\.pcmt\.'
patt += r'[A-D]{1,4}\.[X-Y]{1,2}\.dat$'

pcmt = []
for fn in pcmt_raw:
    pcmt += re.findall(patt, fn)


np.savetxt('pcmt.txt', pcmt, fmt='%s')
