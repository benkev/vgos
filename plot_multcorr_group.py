'''
plot_multcorr_group.py

Generate plots of all the band-pols for one or all the stations,
computing the multiple correlation coefficients of each band-pol with respect
to other band-pols to reject the plots with the mult-corr below the threshold.

Arguments (optional):
  -s <a station letter>, like E (or in lower case, e);
  -o <output directory name> 

If called with no arguments, it processes all the stations E, G, H, I, V, and Y.

'''

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os, sys, glob, copy
import itertools as itr
import re
import datetime, time, calendar
from functools import reduce
import getopt


#
# Process command line options
#
stations = ('E', 'G', 'H', 'I', 'V', 'Y') 
plotdir = ''
optlist = getopt.getopt(sys.argv[1:], 's:o:')[0]
for opt, val in optlist:
    if opt == '-s':
        stations = (val.upper(), )  # A letter like E or e, G or g etc.
    elif opt == '-o':
        plotdir = val
        if plotdir[-1] != '/':
            plotdir += '/'
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)


bp_sym = ['AX', 'AY', 'BX', 'BY', 'CX', 'CY', 'DX', 'DY']
nbp = len(bp_sym)
bp_patt = re.compile(r'[A-D][X-Y]')  # Regex pattern to catch bandpols, like AX

#fields = [('expcode', 'S4'),  ('expname', 'S6'),
#          ('a', 'S2'), ('b', 'S2'), ('c', 'S2'),  ('d', 'S2')]

#
# Dtype for the lines of data in files
# /data/geodesy/CCCC/pcc_datfiles[_jb]/bandmodel.EENNNN.S.C.P.dat
# on the demi.haystack.mit.edu server.
#
fields = [('year', 'i2'),  ('doy', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('phase_midband', 'f8'), ('phase_dc', 'f8'),
          ('delay_ps', 'f8'), ('phase_rmse', 'f8'),
          ('scan', 'S16'), ('source', 'S16'), ('station', 'S2'),
          ('azimuth', 'f8'), ('elevation', 'f8')]

#
# Create hotr, an inverted 'hot' colormap with frduced dynamic range
# (From white throuhg yellow to dense red)
#
cmhot_r = plt.cm.get_cmap('hot_r')
hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')
# hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')
cmYlGn = plt.cm.get_cmap('YlGn')
rYlGn = ListedColormap(cmYlGn(np.linspace(0.0, 0.6, 256)), name='rYlGn')



# ========================================================================
# Read file stationX.txt in the dline list. Remove the blank lines.
# Put into bad_bandpols lists of those  band-polarizations not listed in 
# the stationX.txt file.
# Put into good_bandpols lists of those  band-polarizations listed in 
# the stationX.txt file.
#

for st in stations:
    # fout = open('multcorr_station'+ st +'.txt', 'w')

    with open('station'+ st +'.txt') as fh: 
        dlines = fh.readlines()

    #sys.exit(1)

    #
    # Remove empty lines
    #
    while '\n' in dlines: # Remove empty lines
        dlines.remove('\n')

    print('Station ' + st)

    nexperm = len(dlines)
    bad_bandpols =  [] # List of lists of BAD band-polarizations
    good_bandpols = [] # List of lists of GOOD band-polarizations
    bad = []          # List of lines with the missing data

    for iexperm in range(nexperm):
        bad_bp = ['AX', 'AY', 'BX', 'BY', 'CX', 'CY', 'DX', 'DY'] # Full set
        dl = dlines[iexperm]

        if ('<' in dl) or ('?' in dl) or ('>' in dl):
            bad.append(iexperm)
            #continue # ================================================== >>>

        # if ('-' in dl) or ('!' in dl):
        #     bad_bandpols.append(bad_bp)  # All band-pols are bad
        #     #continue # ================================================== >>>

        #dl = dl[12:-1]        # Leave only chan-pol data
        #bplist = dl.split(',')

        bplist = bp_patt.findall(dl)  # Get a list of the bandpols
        bplist.sort()
        # Leave in bad_bp only those band-pols absent in the dl line 
        for bp in bplist:
            bad_bp.remove(bp)
        
        bad_bandpols.append(bad_bp)
        good_bandpols.append(bplist)

   

    #sys.exit(1)

    # for idx in bad[::-1]:
    #     dlines.pop(idx)

    #sys.exit(1)

    #========================================================================
    

    for idl in range(nexperm):                       # dl in dlines:
        dl = dlines[idl]
        exc = dl[:4]          # Experiment code       
        exn = dl[5:]          # Experiment name (with the end of line)

        #
        # If exn starts with experiment name like VT9133, leave only it in exn
        # Otherwise exn contains something like '<Did not run>'
        #
        vtnnnn = re.findall('^[A-y]{2}[0-9]{4}', exn)
        if vtnnnn:
            exn = vtnnnn[0]
        else:
            exn = exn[:-1]

        wrline1 = exc + ' ' + exn + '  '

        for ibp in range(nbp):         # nbp = 8 bandpols
            bp = bp_sym[ibp]
            if bp in good_bandpols[idl]:
                wrline1 += '    ' + bp + '    '
            else:
                wrline1 += '          '

        # fout.write(wrline1 + '\n')

        ddir = '/data/geodesy/' + exc + '/pcc_datfiles'
        if os.path.isdir(ddir):
            datdir = copy.copy(ddir)
        elif os.path.isdir(ddir + '_jb'):
            datdir = ddir + '_jb'
        else:
            wrline2 = 'Neither "'+ddir+'" nor "'+ddir+ \
                     '_jb" exists. \n     Skipped.'
            print(wrline2)
            # fout.write('++++ '+wrline2 + '\n')
            continue # ============================================== >>>

        fnames = glob.glob(datdir + '/bandmodel.??????.' + st + '.?.?.dat')

        if len(fnames) == 0: # No data for the station st
            continue # ============================================== >>>

        fnames.sort()

        nfile = len(fnames)

        figname = plotdir + 'bandpol_st_' + st + '_' + exc + '_' + exn + \
                  '_ABCD_XY.png'

        #
        # Read all the 8 channel data into the array list datlist.
        # In case the datasets have different lengths, find the minimum length
        # min_ndat to afterwards truncate the dataset arrays in datlist 
        #
        min_ndat = 1000000000
        datlist = []
        ndats = np.zeros(nfile, dtype=int)
        for ix in range(nfile):
            dat_ix = np.loadtxt(fnames[ix], dtype=fields)
            ndats[ix] = dat_ix.shape[0]
            if dat_ix.shape[0] < min_ndat:
                min_ndat = dat_ix.shape[0]
            datlist.append(dat_ix)

        for ix in range(nfile):  # Reduce the arrays in datlist to minimum size
            if ndats[ix] != min_ndat:
                datlist[ix] = datlist[ix][:min_ndat]

        dat = np.array(datlist)  # Convert list of struct. arrays to array
        ndat = np.size(dat,1)    # Height of a column
        delps = np.zeros((nfile,ndat), dtype=float)  # Cable delays (ps)

        for ix in range(nfile):
            delps[ix,:] = dat[ix]['delay_ps']

        #
        # Assume time data the same for all bands. Use time from AX (ie 0-th).
        #
        tyear =   dat[0]['year'].astype(int)
        tdoy =    dat[0]['doy'].astype(int)
        thour =   dat[0]['hour'].astype(int)
        tminute = dat[0]['minute'].astype(int)
        tsecond = dat[0]['second'].astype(int)

        datim0 = datetime.datetime(int(tyear[0]), 1, 1) + \
                 datetime.timedelta(int(tdoy[0]) - 1)
        ttuple0 = datim0.timetuple()
        tstamp0 = time.mktime(ttuple0)

        exp_doy_time = datim0.strftime('%j, %H:%M:%S')


        t_sec = np.asarray(tsecond, dtype=float)

        for itim in range(ndat):
            datim = datetime.datetime(int(tyear[itim]), 1, 1, \
                                      int(thour[itim]), int(tminute[itim]), \
                                      int(tsecond[itim])) + \
                    datetime.timedelta(int(tdoy[itim]) - 1)
            ttuple = datim.timetuple()
            tstamp = time.mktime(ttuple)
            t_sec[itim] = tstamp

        t_hr = (t_sec - tstamp0)/3600.    # Time in hours  






        nbandpol = nfile # Assume there are 8 files for each band and pol

        #
        # Compute 8x8 correlation matrix of rows of delps[8,878]
        #
        Rxx_full = np.corrcoef(delps)

        #
        # Compute the multiple correlation coefficients for
        # every bandpol
        #          
        R_mult2 = np.zeros(nbandpol, dtype=float)  # Squared mult corrs
        R_mult = np.zeros(nbandpol, dtype=float)   # Mult corr coefficients
        bprange = np.arange(nbandpol) 
        bpxy = np.array([bp for bp in itr.combinations(bprange,2)], dtype=int)
        bpx = bpxy[:,0]
        bpy = bpxy[:,1]
        for ibp in range(nbandpol):
            #
            # The ibp-th bandpol is assumed an independent variable.
            # Rxx is the correlation matrix of 7 other bandpols;
            # cor is the vector of cross-correlations of each of the 7 bandpols
            # with the ibp-th bandpol.
            #
            Rxx = np.delete(np.delete(Rxx_full, ibp, axis=0), ibp, axis=1)
            cor = np.delete(Rxx_full[ibp,:], ibp, axis=0)
            invRxx = la.inv(Rxx)
            # R2 = multi_dot([cor, Rxx, cor])
            # R2 = cor.dot(Rxx).dot(cor)
            R_mult2[ibp] = reduce(np.dot, [cor, invRxx, cor])
        R_mult = np.sqrt(R_mult2)
        R_percent = 100.*R_mult
        
        #
        # Write station file with multcorr values
        #
        wrline2 = 5*' ' + exn + ' '
        for ibp in range(nbandpol):
            wrline2 += ' {:7.4f}  '.format(100*R_mult[ibp])

        # fout.write(wrline2 + '\n')



        #
        # Create a plot for each band/pol we have data for (on a 2X4 grid)
        #
        fig = plt.figure(figsize=(8.5,11))
        fig.suptitle("Exp. " + exn + " (code " + exc + \
                     "), Station " + st + \
                     ". Delay trend for bands ABCD:XY and R_multcorr")
        strpol =  'XY'
        strband = 'ABCD'

        iplot = 0
        for ibp in range(nbandpol):
            ip = ibp % 2        # 0 1 0 1 0 1 0 1
            ib = ibp // 2        # 0 0 1 1 2 2 3 3
            iplot = iplot + 1   # subplot index starts from 1

            band_pol = strband[ib] + ':' + strpol[ip]

            ax = plt.subplot(4, 2, iplot)
            ax.plot(t_hr, delps[ibp,:], 'b.', markersize=3, clip_on=False)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            x_marker = xmin + 0.93*(xmax - xmin)
            y_marker = ymin + 0.9*(ymax - ymin)
            x_text = xmin + 0.73*(xmax - xmin)
            y_text = ymin + 0.87*(ymax - ymin)
            icol = int(R_mult[ibp]*255)
 
            ax.plot(x_marker, y_marker, marker='s', markersize=20, \
                     markerfacecolor=rYlGn.colors[icol])
            ax.text(x_text, y_text, '%5.2f' % R_percent[ibp])

            if bp_sym[ibp] in bad_bandpols[idl]:
                x_text2 = xmin + 0.05*(xmax - xmin)
                y_text2 = ymin + 0.87*(ymax - ymin)
                ax.text(x_text2, y_text2, 'Rejected', color='r')

            ax.set_ylabel(band_pol + " delay (ps)")
            ax.grid(True)

            if iplot == 7 or iplot == 8: # bottom two plots
                ax.set_xlabel("hours since doy " + exp_doy_time)

        fig.tight_layout()
        
        fig.subplots_adjust(top=0.94)
        fig.savefig(figname)

        #plt.show()
        plt.close(fig)
        
        #raise SystemExit


    # fout.close()



