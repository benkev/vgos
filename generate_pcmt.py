'''
generate_pcmt.py

This program generates pcmt files using pcc_select.py with the bandpol
selections generated by select_bandpols.py and saved in the file 
'selections_st_<station(s)>_m<median-threshold>.txt'.

Arguments (positional):
    <band-pol-selections-file-name>
    <directory-to-save-pcmt>

Example:
  %run generate_pcmt.py pltE_m0.97/selections_st_E_m0.97.txt pcmt_E_m0.97/

'''

import numpy as np
import re, os, sys
import getopt

sfname = sys.argv[1]   # File name with the bandpol selections
odname = sys.argv[2]   # Output directory name to store pcmt files

#
# Process command line options
#
# if sys.argv[1:] == []: # Print help text and exit if no command line options
#     print(help_text)
#     raise SystemExit

# optlist = getopt.getopt(sys.argv[1:], 'm:s:d:o:phxa')[0]

# for opt, val in optlist:
#     if opt == '-h':  # Print help text and exit if there is '-h' among options
#         print(help_text)
#         raise SystemExit
# 

# pcmt = np.loadtxt('pcmt.txt', dtype='a128', delimiter='\n')

with open('pcmt.txt') as f: pcmt = f.readlines()

#
# Get the codes for experiments with existing pcmt files
#
# To not include the slashes, use "matching group". 
# ? means zero or one occurrence.
#
patt = r'\/(\d{4}?)\/'   
exp_codes = [re.findall(patt, pdir)[0] for pdir in pcmt]
exp_codes = map(int, exp_codes)
exp_codes = np.unique(exp_codes)
exp_codes.sort() # Experiment codes for the existing pcmt files

sfdir = os.path.dirname(sfname)

sel_lines = np.loadtxt(sfname, dtype='a256', delimiter='\n')

#for sl in sel_lines[:100]:

#
# For the experiments with existing pcmt files only,
# create pcmt files based on select_bandpols.py selections
#
for exc in exp_codes:

    sel_found = False
    for isl in range(len(sel_lines)):
        sl = sel_lines[isl]
        if sl[:4] == str(exc):  # and sl[12:14] != 'jb':
            sel_found = True
            break

    if not sel_found:
        continue  # ======================================================= >>>
    
    sels = re.findall('[EGHIVY]:(?:[A-D][XY],*){,8}', sl)

    # sels: a list like 
    # ['E:AX,AY,BX,BY,CX,CY,DX,DY',
    #  'G:BX,BY,CX,CY,DX,DY',
    #  'H:AX,AY,BX,BY,CX,CY,DX,DY',
    #  'I:CX,CY,DX,DY',
    #  'V:AX,AY,BX,BY,CX,CY,DX',
    #  'Y:AX,AY,CX,CY']

    sels_line = ' '.join(sels)      # Put the list into a blank-separated line
    sels_line += ' '
    excode = sl[:4]
    exname = sl[4:12]
    pcc_select_cmd =  'pcc_select.py -e' + exname
    pcc_select_cmd += '-d /data/geodesy/' + excode + '/pcc_datfiles/ '
    pcc_select_cmd += '-s ' + sels_line + ' '
    pcc_select_cmd += '-o ' + odname

    print(pcc_select_cmd)
    
    os.system(pcc_select_cmd)


# raise SystemExit

