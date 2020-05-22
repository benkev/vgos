'''
plot_multcorr_all.py
'''
help_text = '''
Generates plots of the band-pols for one experiment on one station,
Computs the multiple correlation coefficients of each band-pol with respect
to other band-pols.
Computes medians for the rows (or columns) of the correlation matrix.
The plots with the mult-corr below its threshold (-t or 90.) or with the 
correlation median below its threshold (-m or 0.5) are marked as "Rejected". 

Arguments:
  -t <threshold>              for multiple correlation coefficient, 0. to 100.
  -m <threshold>              for correlation median, -1 to 1.
  -s <a station letter>, like E, G, H ... (or in lower case, e, g, h ...);
  -d <pcc_datfiles directory>       like /data/geodesy/3686/pcc_datfiles
  -o <output directory name>        where .png graphs and .txt logs are saved
  -x                          plot x-correlation matrix
  -p                          show plot in X-window
  -a                          Make .png plots and .txt files for all available
                              data under directory in -d (like -d /data/geodesy)
                              If one or more stations are given in -s, like
                              -s E or -s VIGH, only data for those stations are
                              plotted and saved in -o directory.
                              If -a is present, -p is ignored (for too many 
                              wind would be open).
  -h                          print this text'''


import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from phasecal import mult_corr
import os, sys, glob, copy
import itertools as itr
import re
import datetime, time, calendar
from functools import reduce
import getopt


threshold_mulcor = 90.
threshold_median = 0.5

stations = ('E', 'G', 'H', 'I', 'V', 'Y') 
bp_sym = ['AX', 'AY', 'BX', 'BY', 'CX', 'CY', 'DX', 'DY']
nbp = len(bp_sym)
bp_patt = re.compile(r'[A-D][X-Y]')  # Regex pattern to catch bandpols, like AX

fields = [('year', 'i2'),  ('doy', 'i2'),
          ('hour', 'i2'), ('minute', 'i2'), ('second', 'i2'),
          ('phase_midband', 'f8'), ('phase_dc', 'f8'),
          ('delay_ps', 'f8'), ('phase_rmse', 'f8'),
          ('scan', 'S16'), ('source', 'S16'), ('station', 'S2'),
          ('azimuth', 'f8'), ('elevation', 'f8')]


#
# Process command line options
#
if sys.argv[1:] == []: # Print help text and exit if no command line options
    print(help_text)
    raise SystemExit

optlist = getopt.getopt(sys.argv[1:], 't:s:d:o:phxa')[0]

raise SystemExit

for opt, val in optlist:
    if opt == '-h':  # Print help text and exit if there is '-h' among options
        print(help_text)
        raise SystemExit

datadir = ''
station = ''
outdir = ''
nbandpol = 0
plot_graph = False
datadir_exists = True
station_data_exist = True
plot_xcorrmx = False

for opt, val in optlist:
    if opt == '-t':
        t_0 = float(val)
        if t_0 < 0.:
            threshold_mulcor = 0.
        elif t_0 > 100.:
            threshold_mulcor = 100.
        else:
            threshold_mulcor = t_0
    if opt == '-m':
        m_0 = float(val)
        if m_0 < -1.:
            threshold_median = -1.
        elif m_0 > 1.:
            threshold_median = 1.
        else:
            threshold_mulcor = m_0
    if opt == '-s':
        station = val.upper()     # A letter like E or e, G or g etc.
    elif opt == '-d':
        datadir = copy.copy(val)
    elif opt == '-o':
        outdir = copy.copy(val)   # Like /data/geodesy/3671/pcc_datfiles_jb
        if outdir[-1] != '/':
            outdir += '/'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    elif opt == '-x':
        plot_xcorrmx = True
    elif opt == '-p':
        plot_graph = True

if station == '':
    print('ERROR: Station must be specified in option -s.')
    station_data_exist = False
elif not station in stations:
    print('ERROR: Station ' + station + ' not in ' + str(stations))
    station_data_exist = False

if datadir == '':
    print('ERROR: Data directory must be specified in option -d.')
    station_data_exist = False
elif os.path.isdir(datadir):
    fnames = glob.glob(datadir + '/bandmodel.??????.' + station + '.?.?.dat')
    nbandpol = len(fnames)
    if nbandpol == 0:          # No data for the station 
        print('ERROR: No data for station ' + station + ' on path ' + datadir)
        station_data_exist = False
    else:
        fnames.sort()
else:
    print('ERROR: Path ' + datadir + ' does not exist.')
    datadir_exists = False


if '' in [datadir, station]:
    raise SystemExit
if False in [datadir_exists,station_data_exist]: 
    raise SystemExit

#
# Extract from path the experiment code like /3686/
#
exc = re.findall('\/[0-9]{4}\/', datadir)
exc = exc[0][1:-1]

#
# Extract from data file name the experiment name like .VT9050.
#
exn = re.findall('\.[A-y]{2}[0-9]{4}\.', fnames[0])
exn = exn[0][1:-1]

outname = outdir + 'bandpol_st_' + station + '_' + exc + '_' + exn + \
          '_ABCD_XY'
figname = outname + '.png'
txtname = outname + '.txt'

outname2 = outdir + 'xcorr_matrix_st_' + station + '_' + exc + '_' + exn + \
          '_ABCD_XY'
figname2 = outname2 + '.png'
# txtname2 = outname2 + '.txt'


#
# Create hotr, an inverted 'hot' colormap with frduced dynamic range
# (From white throuhg yellow to dense red)
#
# cmhot_r = plt.cm.get_cmap('hot_r')
# hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')

#
# Create rYlGn, a YlGn (Yellow-Green) colormap with frduced dynamic range
#
cmYlGn = plt.cm.get_cmap('YlGn')
rYlGn = ListedColormap(cmYlGn(np.linspace(0.0, 0.6, 256)), name='rYlGn')


#
# Read all the 8 channel data into the array list datlist.
# In case the datasets have different lengths, find the minimum length
# min_ndat to afterwards truncate the dataset arrays in datlist 
#
min_ndat = 1000000000
datlist = []
ndats = np.zeros(nbandpol, dtype=int)
for ix in range(nbandpol):
    dat_ix = np.loadtxt(fnames[ix], dtype=fields)
    ndats[ix] = dat_ix.shape[0]
    if dat_ix.shape[0] < min_ndat:
        min_ndat = dat_ix.shape[0]
    datlist.append(dat_ix)

for ix in range(nbandpol):  # Reduce the arrays in datlist to minimum size
    if ndats[ix] != min_ndat:
        datlist[ix] = datlist[ix][:min_ndat]

dat = np.array(datlist)  # Convert list of struct. arrays to array
ndat = np.size(dat,1)    # Height of a column
delps = np.zeros((nbandpol,ndat), dtype=float)  # Cable delays (ps)

#
# The cable delays for all band-pols are put in 2D appay delps
#
for ix in range(nbandpol):
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


#
# Compute nbandpol X nbandpol correlation matrix of rows of delps[nbandpol,:]
#
Rxx_full = np.corrcoef(delps)

#
# Compute medians for the rows (or columns) of the correlation matrix
# Low median points at too weak correlation of a selected bandpol with
# other bandpols.
#
Rxx_full_0 = np.copy(Rxx_full)
for i in range(nbandpol):      # Replace ones at the diagonal with zeros
    Rxx_full_0[i,i] = 0.

corr_median = np.median(Rxx_full, axis=0)


#
# Compute the multiple correlation coefficients for every bandpol
#          
R_mult = mult_corr(Rxx_full)
R_percent = 100.*R_mult



#
# Write log file with multcorr values
#
fout = open(txtname, 'w')

wrl1_1 = '#\n# Mulpiple Correlation Coefficients (%). Station ' + station + \
         ', Exp. ' + exn + ', Code ' + exc + '\n#\n'
wrl1_2 = '#'

fout.write(wrl1_1 + '\n')

wrline1 = '# ' + exc + ' ' + exn + ' '
for ibp in range(nbp):         # nbp = 8 bandpols
    wrline1 += '    ' + bp_sym[ibp] + '  '

wrline2 = 15*' '
for ibp in range(nbandpol):
    wrline2 += ' {:5.2f}  '.format(R_percent[ibp])

fout.write(wrline1 + '\n')
fout.write(wrline2 + '\n\n')

#
# Save the cross-correlation matrix in file
#
wrl2_1 = '#\n# Cross-Correlation Matrix. Station ' + station + \
         ', Exp. ' + exn + ', Code ' + exc + '\n#\n'
wrl2_2 = 14*' '
for ibp in range(nbp):         # nbp = 8 bandpols
    wrl2_2 += '    ' + bp_sym[ibp] + '  '

fout.write(wrl2_1)
fout.write(wrl2_2 + '\n')

for iy in range(nbandpol):
    wrl = 11*' ' + bp_sym[iy] + ' '
    for ix in range(nbandpol):
        wrl += ' {:6.3f} '.format(Rxx_full[iy,ix])
    fout.write(wrl + '\n')
fout.write('\n')

#
# Save the cross-correlation medians in file
#
wrl = '      Medians '
for ix in range(nbandpol):
    wrl += ' {:6.3f} '.format(corr_median[ix])

fout.write(wrl + '\n\n')

fout.close()



#
# Create a plot for each band/pol we have data for (on a 2X4 grid)
#
fig = plt.figure(figsize=(8.5,11))
fig.suptitle("Exp. " + exn + " (code " + exc + \
             "), Station " + station + \
             ". Delay trend for bands ABCD:XY and R_multcorr")
strpol =  'XY'
strband = 'ABCD'

iplot = 0
for ibp in range(nbandpol):
    ip = ibp % 2         # 0 1 0 1 0 1 0 1
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
    
    if (R_percent[ibp] < threshold_mulcor) or \
       (corr_median[ibp] < threshold_median):
        x_text2 = xmin + 0.05*(xmax - xmin)
        #x_text3 = xmin + 0.30*(xmax - xmin)
        #y_text2 = ymin + 0.87*(ymax - ymin)
        ax.text(x_text2, y_text, 'Rejected', color='r')

    # Print correlation median
    x_text3 = xmin + 0.30*(xmax - xmin)
    ax.text(x_text3, y_text, '%5.2f' % corr_median[ibp], color='k')


    ax.set_ylabel(band_pol + " delay (ps)")
    ax.grid(True)

    if iplot == 7 or iplot == 8: # two bottom plots
        ax.set_xlabel("hours since doy " + exp_doy_time)

fig.tight_layout()

fig.subplots_adjust(top=0.94)
fig.savefig(figname)

if plot_xcorrmx:
    #
    # Do not plot the half of cross-correlation matrix under the diagonal
    #
    Rxx_nan = np.copy(Rxx_full)
    for ix in range(nbandpol):
        Rxx_nan[ix,:(ix+1)] = np.NaN

    #
    # Use inverted 'hot' colormap with reduced dynamic range
    # (From white throuhg yellow to dense red)
    #
    cmhot_r = plt.cm.get_cmap('hot_r')
    hotr = ListedColormap(cmhot_r(np.linspace(0.0, 0.7, 256)), name='hotr')
    n0_7 = np.arange(8)  # Tick values for 

    fig2 = plt.figure(figsize=(6,5));
    ax2 = plt.subplot(111)

    xcorimg = ax2.imshow(Rxx_nan, interpolation='none', cmap=hotr);
    #xcorimg = ax2.imshow(Rxx_nan, interpolation='none', cmap=rYlGn);
    #plt.pcolormesh(Rxx_full, cmap=plt.cm.jet, offset_position='data');
    ax2.set_xticks(n0_7)
    ax2.set_xticklabels(bp_sym)
    ax2.tick_params(axis='x', labeltop='on')
    ax2.set_yticks(n0_7)
    ax2.set_yticklabels(bp_sym)
    #ax2.grid(1)
    fig2.colorbar(xcorimg, shrink=0.8)
    fig2.text(0.1, 0.95, 'Cross-Correlation Matrix. Station ' + station +
                 ', Exp. ' + exn + ', Code ' + exc)

    #fig.tight_layout()
    #fig.subplots_adjust(top=0.94)

    fig2.savefig(figname2)


if plot_graph:
    fig.show()
    fig2.show()
else:
    plt.close(fig)
    plt.close(fig2)
        

# fout.close()



