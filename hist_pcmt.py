'''
hist_pcmt.py

Plots 3x4=12 histograms of the rms difference between PCMT files.
Each plot has 3x4=12 histograms for a station at the same median threshold 
value. The histograms are built for all the stations available.

Also, it plots graphs of average RMS of difference (ps) vs median thresholds.

The rms data are read from files pcmt_stat.txt in directories pcmt_m0.xx/

The first parameter on the command line is the output directory name.

The set of 12 thresholds is specified on the command line. The threshold values
can be either numbers in [0..1] at the precision of two decimals after the dot
or just 2-digit percent values. I.e., the following numbers are considered
equivalent:
  0.78 or .78 or 78

Example:

python hist_pcmt.py hist_pcmt1 50 70 90 91 92 93 94 95 96 97 98 99

or, which is the same:

python hist_pcmt.py hists/ .5 0.7 90 0.91 .92 .93 .94 0.95 96 97 98 99

'''

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import re, os, sys, glob, copy
import getopt

# fields = [('year', 'i4'),  ('month', 'i2'),  ('day', 'i2'),
#           ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
#           ('delay_s', 'f8'), ('source', 'a9'), ('delay_ps', 'a9')]

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

st1name = {
    'E' : 'Westford',
    'F' : 'Westford',
    'G' : 'Goddard',
    'H' : 'Kokee',
    'V' : 'Wettzell',
    'Y' : 'Yebes',
    'I' : 'Ishioka',
    'S' : 'Onsala-NE',
    'T' : 'Onsala-SW',
    'M' : 'MacDonald'
}

st2name = {
    'wf' : 'Westford',
    'w2' : 'Westford',
    'gs' : 'Goddard',
    'k2' : 'Kokee',
    'ws' : 'Wettzell',
    'yj' : 'Yebes',
    'is' : 'Ishioka',
    'oe' : 'Onsala-NE',
    'ow' : 'Onsala-SW',
    'mg' : 'MacDonald'
}


outdir = sys.argv[1]
if outdir[-1] != '/': outdir += '/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)

thrs = sys.argv[2:]


if len(thrs) == 0:
    print("Usage:\npython hist_pcmt.py <thr1> <thr2> <thr3> ... <thrN>")
    raise SystemExit

#raise SystemExit

#
# Reduce the cmd line numbers to the proper format like 0.xx
# and save them in thresholds list
#
thresh = []
str_thresh = []
nth = 0
for th in thrs:
    thf = float(th)
    if thf > 1.: thf = 0.01*thf
    str_thresh.append('%4.2f' % thf)
    thresh.append(thf)
    nth = nth + 1
#    if nth >= 12: # Leave only maximum 12 thresholds 
#        break     # ====================================================== >>>

#
# Dictionary to store rms-es for individual stations for individual thresholds
#
rmsd_allst = { 
    'wf' : OrderedDict(), # Westford
    'w2' : OrderedDict(), # Westford2 (defunct)
    'gs' : OrderedDict(), # GGAO12M (Goddard)
    'k2' : OrderedDict(), # Kokee
    'ws' : OrderedDict(), # Wettzell
    'yj' : OrderedDict(), # Yebes
    'is' : OrderedDict(), # Ishioka
    'oe' : OrderedDict(), # Onsala-Northeast
    'ow' : OrderedDict(), # Onsala-Southwest
    'mg' : OrderedDict()  # MGO (MacDonald)
}

for th in str_thresh:
    #
    # Find and read in the relevant file
    #
    indir = 'pcmt_m' + th + '/'
    if not os.path.isdir(indir):
        print('WARNING: no such directory: ' + indir + '. Skipping.')
        continue   # ====================================================== >>>
    infile = indir + 'pcmt_stat.txt'
    if not os.path.isfile(infile):
        print('WARNING: no such file: ' + infile + '. Skipping.')
        continue   # ====================================================== >>>
    
    with open(indir + 'pcmt_stat.txt') as fin: 
        exp_st_rms = fin.readlines() # Lines with experiment, station, and rms
    exp_st_rms = exp_st_rms[:-1] # Throw away last line with 'mean rms:'
        
    #
    # Loop over the lines and distribute the rms-es to the stations
    #
    for esr in exp_st_rms:
        st = esr[6:8]       # Station 2-char code
        rms = float(esr[8:-1])
        if th not in rmsd_allst[st]:
            rmsd_allst[st][th] = [rms]  # Create list with the first rms
        else:
            rmsd_allst[st][th].append(rms)

#
# Remove unused stations from rmsd_allst[st][th] and put the data into 
# the rmsd[st][th] dict converting lists into np.appays.
#
rmsd = OrderedDict()
for st in rmsd_allst.keys():
    if rmsd_allst[st] != OrderedDict():
        rmsd[st] = OrderedDict()

for st in rmsd.keys():
    for th in rmsd_allst[st].keys():
        rmsd[st][th] = copy.copy(rmsd_allst[st][th])


#raise SystemExit

#
# RMS Averages:
#
rmsd_avg = OrderedDict()

#raise SystemExit

# outdir = 'hist_pcmt/'
# if not os.path.isdir(outdir):
#    os.mkdir(outdir)

#
# Plot average RMS of difference (ps) vs median threshold 
#
fig1 = plt.figure(figsize=(10,8))
fig1.suptitle('Average RMS of Difference (ps) vs Median Threshold', fontsize=18)
# fig1.suptitle('Station "' + st + '". RMS of Diff (ps) ' \
#                  'for 12 Median Threshs., %d Values' % nrms, fontsize=18)

iavrms = 0

for st in rmsd.keys():
    #
    # Assume the number of median values is the same for all thresholds
    #
    nrms = len(rmsd[st][str_thresh[0]])

    #
    # Find maximal rms for the station st to set common histogram x-scale.
    # Remove too large rms values.
    #
    st_rms_max = -1e6
    for th in rmsd[st].keys():
        rmsl = rmsd[st][th]
        rms_max1 = max(rmsl)
        while rms_max1 > 50:         # Remove spurious values
            rmsl.remove(rms_max1)
            rms_max1 = max(rmsl)
        if st_rms_max < rms_max1:
            st_rms_max = rms_max1

    #
    # For the station st, put the rms averages for each threshold level 
    # into rmsd_avg[st].
    # Find minimum among the rms averages avg_rms_min 
    # and its threshold key th_min.
    #
    rmsd_avg[st] = OrderedDict()
    avg_rms_min = 1e6 
    for th in rmsd[st].keys():
        avg_rms = np.mean(rmsd[st][th])
        rmsd_avg[st][th] = avg_rms
        if avg_rms < avg_rms_min:
            avg_rms_min = avg_rms
            th_min = th
            fth_min = float(th_min)

    print('\nStation "%s"; st_rms_max = %f' % (st, st_rms_max))

    if nrms <= 2: # Do not plot histogram for too few values
        continue  # ================================================== >>>

    #
    # Plot average RMS of difference (ps) vs median thresholds
    #
    avrms = [val for val in rmsd_avg[st].values()]

    ax = fig1.add_subplot(2,3,iavrms+1)
    ax.plot(thresh, avrms, lw=2, color='blue')
    ax.plot(thresh, avrms, 'bo', markersize=4, color='black')
    ax.plot([fth_min], [avg_rms_min], 'r*', markersize=10)
    ax.text(0.25, 0.9, 'Station "' + st + '"', transform = ax.transAxes, \
            fontsize=18)
    ax.text(0.15, 0.3, 'Optimum at %s' % th_min, fontsize=14, \
            transform = ax.transAxes)
    ax.grid(1)
    fig1.tight_layout(rect=[0, 0, 1, 0.96])
    iavrms = iavrms + 1

    #
    # Plot one 3x4-figure per station 
    #
    # Plot 4x3=12 histograms of the RMS of difference (ps) between
    # two pcmt files for 12 Median Threshs
    #
    fig2 = plt.figure(figsize=(10,12))
    fig2.suptitle('Station "' + st + '". RMS of Difference (ps) ' \
                  'for 12 Median Threshs., %d Values' % nrms, fontsize=14)
    ihist = 0

    for th in rmsd[st].keys():
        rms = np.array(rmsd[st][th], dtype=float)
        ax = fig2.add_subplot(4,3,ihist+1)
        ax.hist(rms, range=[0., np.ceil(st_rms_max)], alpha=0.7)
        ax.grid(1)
        txtcol = 'red' if th == th_min else 'black' 
        ax.text(0.5, 0.9, 'avg. rms  ' + '%5.2f' % rmsd_avg[st][th], \
                transform = ax.transAxes, color=txtcol)
        ax.text(0.5, 0.8, 'm. thresh. ' + th, transform = ax.transAxes, \
                color=txtcol)
        fig2.tight_layout(rect=[0, 0, 1, 0.96])

        ihist = ihist + 1

        if ihist >= 12:
            break   # ==================================================== >>>


    # fig1.show()
    # fig2.show()
    # raise SystemExit

    fig2.savefig(outdir + 'hists_' + st + '.png')
    plt.close(fig2)

fig1.text(0.69, 0.25, 'X axis: median threshold', fontsize=16)
fig1.text(0.69, 0.20, 'Y axis: average rms', fontsize=16)

fig1.savefig(outdir + 'avg_rms_vs_thresh.png')
#fig1.show()
plt.close(fig1)

