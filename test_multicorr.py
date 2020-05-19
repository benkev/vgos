'''
test_multicorr.py
'''
help_text = '''
Test the multiple correlation of correlated and anti-correlated random
sequences

Arguments:
  -o <output directory name>        where .png graphs and .txt logs are saved
  -p                                show plot (along with saving in .png file)
  -h                                print this text
'''

import numpy as np
from numpy.random import normal
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from functools import reduce
import os, sys, glob, copy
import getopt


def mult_corr(Rxx_full, bp_good=None):
    '''
    Compute multiple correlation coefficients for selected
    bandpols, whose indices are in the sequence bp_good.

    For each i-th bandpol the square of multiple correlation coefficient 
    is computed  as vector-matrix-vector product
      R^2_mult[i] = C.T * Rxx^-1 * C,
    where 
      C = cor, the vector of cross-correlations of each of the bandpols
               except the i-th
      C.T is C transposed,
      Rxx is the matrix of cross-correlations between the band-pol series
             with the i-th bandpol excluded.
          
    Input parameters:
      Rxx_full[nbandpol,nbandpol] correlation matrix of the cable delay 
                 time sequences in the rows of delps[nbandpol,:]
      bp_good: list of indices of good band-pols

    Returns:
      R_mult, array of multiple correlation coefficients
    '''
    if bp_good is None:
        Rxx_good = np.copy(Rxx_full)
    else:
        #
        # From matrix Rxx_full, extract a sub-matrix Rxx_good
        # having only columns and rows with the indices from list bp_good
        #
        cols = list(bp_good)          # Like [ 2,   3,   4,   5,   7 ]
        rows = [[i] for i in cols]    # Like [[2], [3], [4], [5], [7]]

        Rxx_good = Rxx_full[rows,cols]


    nbandpol_good = np.size(Rxx_good, 0)
    R_mult2 = np.zeros(nbandpol_good, dtype=float)  # Squared mult corrs
    R_mult = np.zeros_like(R_mult2)                 # Mult corr coefficients

    for ibp in range(nbandpol_good):
        #
        # The ibp-th bandpol is assumed an independent variable.
        #
        # Rxx is the correlation matrix of the bandpols except the ibp-th one
        # obtained by scratching off ibp-th row and ibp-th column from Rxx_good
        Rxx = np.delete(np.delete(Rxx_good, ibp, axis=0), ibp, axis=1)

        # cor is the vector of cross-correlations of each of the bandpols,
        # except the ibp-th one, with the ibp-th bandpol.
        cor = np.delete(Rxx_good[ibp,:], ibp, axis=0)

        invRxx = la.inv(Rxx)
        R_mult2[ibp] = reduce(np.dot, [cor, invRxx, cor])  # = C.T * Rxx^-1 * C

    R_mult = np.sqrt(R_mult2)
    
    return R_mult

# np.set_printoptions(precision=4)

#
# Process command line options
#
if sys.argv[1:] == []: # Print help text and exit if no command line options
    print(help_text)
    raise SystemExit

optlist = getopt.getopt(sys.argv[1:], 'o:ph')[0]

for opt, val in optlist:
    if opt == '-h':  # Print help text and exit if there is '-h' among options
        print(help_text)
        raise SystemExit


outdir = ''
plot_graph = False

#datadir = ''
#station = ''
#nbandpol = 0
#datadir_exists = True
#station_data_exist = True

for opt, val in optlist:
    # if opt == '-t':
    #     t_0 = float(val)
    #     if t_0 < 0.:
    #         threshold_0 = 0.
    #     elif t_0 > 100.:
    #         threshold_0 = 100.
    #     else:
    #         threshold_0 = t_0
    # if opt == '-s':
    #     station = val.upper()     # A letter like E or e, G or g etc.
    # elif opt == '-d':
    #     datadir = copy.copy(val)
    if opt == '-o':
        outdir = copy.copy(val)
        if outdir[-1] != '/':
            outdir += '/'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    elif opt == '-p':
        plot_graph = True





    
nrprocs = 3
ndat = 1000000
noise = normal(0., 1., (nrprocs,ndat))  # Random processes
rprocs = normal(0., 1., (nrprocs,ndat))  # Random processes


#
# Compute nrprocs*nrprocs correlation matrix of rows of rprocs[nrprocs,ndat]
#
Rxx_full = np.corrcoef(rprocs)

R_mult = mult_corr(Rxx_full)


outname = outdir + 'test_multicorr'
figname = outname + '.png'
txtname = outname + '.txt'

#
# Open log file to save the  multcorr values
#
fout = open(txtname, 'w')

wrline1 = ''
for irpc in range(nrprocs):                         # nbp = 8 bandpols
    wrline1 += '    ' + str(irpc) + '    '
fout.write(wrline1 + '\n')
print(wrline1)
wrline2 = ''
for irpc in range(nrprocs):
    wrline2 += ' {:7.4f}  '.format(R_mult[irpc])
fout.write(wrline2 + '\n')
print(wrline2)

