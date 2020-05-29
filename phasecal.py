'''
phasecal.py

Module providing some useful functions for the phase calibration software.

'''

import numpy as np
import numpy.linalg as la
from functools import reduce


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
      R_mult, array of multiple correlation coefficients, one for each
              bandpol indicated in bp_good.
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



def write_xcorrmx(fout, Rxx_full, station, experiment_number, experiment_code, \
                  bp_sym, ):
    '''
    Save the cross-correlation matrix in file fout.
    '''

    nbandpol = np.size(Rxx_full,0)
    nbp_sym = len(bp_sym)

    wrl2_1 = '#\n# Cross-Correlation Matrix. Station ' + station + \
             ', Exp. ' + experiment_number + ', Code ' + experiment_code + \
             '\n#\n'
    wrl2_2 = 14*' '
    for ibp in range(nbp_sym):         # nbp_sym = 8 bandpols
        wrl2_2 += '    ' + bp_sym[ibp] + '  '

    fout.write(wrl2_1)
    fout.write(wrl2_2 + '\n')

    for iy in range(nbandpol):
        wrl = 11*' ' + bp_sym[iy] + ' '
        for ix in range(nbandpol):
            wrl += ' {:6.3f} '.format(Rxx_full[iy,ix])
        fout.write(wrl + '\n')
    fout.write('\n')

