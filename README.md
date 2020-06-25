# vgos

Unaided selection of good cable delay data.

The coefficient of multiple correlation of each band-pol with respect 
to other seven is used. Also, medians for the rows (or columns) of the
correlation matrix are calcuated 

===============================================================================



select_bandpols.py

help_text = '''
This script selects good band-pols whose rows in the correlation matrix have
medians below the threshold (now 0.5). The rows and columns of bad bandpols are
iteratively removed from the matrix, improving the median values for other
band-pols, until all medians are above the threshold. The correlation matrix,
medians, and multiple correlation coefficients on each iteration (if any) are
saved in text files. 

Unlike plot_multcorr.py, the script can work on all the data (option -a) for
one or several stations (like -s E, -s EGY, or -s IEHV). 

Generates plots of the band-pols for one experiment on one station, or all the
experiments on one or more stations. Computes medians for the rows (or columns)
of the correlation matrix. Computes the multiple correlation coefficients of 
each band-pol with respect to other band-pols.

The plots with with the correlation median below its threshold are marked as 
"Rejected". 

Arguments:
  -m <threshold>              for correlation median, -1 to 1.
  -s <a station letter(s)>, like E, G, H, IEV (or in lower case, e, g, h, iev);
  -d <pcc_datfiles directory>       like /data/geodesy/3686/pcc_datfiles
  -o <output directory name>        where .png graphs and .txt logs are saved
  -x                          plot x-correlation matrix
  -p                          show plot in X-window
  -a                          Make .png plots and .txt files for all available
                              data under directory in -d (like -d /data/geodesy)
                              If more stations are given in -s, like
                              -s EY or -s VIGH, only data for those stations are
                              plotted and saved in -o directory.
                              If -a is present, -p is ignored (for too many 
                              windows would be open).
  -h                          print this text.

Examples:

(1) Select good band-pols in experiment 3658 from station E. Save the band-pol
plot in directory pltE. Save the correlation matrix plot (-x). Show both plots
on the screen (-p). 
The process of successive selection of good band-pols is logged in text file 
bandpol_E_3658_VT8204.txt. 
The plots are saved as bandpol_E_3658_VT8204.png and xcorrmx_E_3658_VT8204.png: 

%run select_bandpols.py -s E -d /data/geodesy/3658/pcc_datfiles_jb/ -o pltE \
                          -p -x

(2) Select good band-pols in all (-a) the experimental data from station V
located under directory /data/geodesy/. 
Save the band-pol plots in directory pltV. Save the correlation matrix plots
(-x). The key -p (show plots on the screen) is ignored when -a is used.  
The processes of successive selection of good band-pols are
logged in text file bandpol_V_<exp.code>_<exp.name>.txt. 
The plots are saved as bandpol_E_<exp.code>_<exp.name>.png and 
xcorrmx_E_<exp.code>_<exp.name>.png.
The diagnostic messages are logged in diagnostics_st_E.txt:

%run select_bandpols.py -s V -d /data/geodesy/ -o pltV -p -x -a

(3) Select good band-pols in all (-a) the experimental data from stations
I and Y. Save the band-pol plots in directory pltV. Save the correlation matrix
plots (-x).  
The processes of successive selection of good band-pols are
logged in text file bandpol_<station>_<exp.code>_<exp.name>.txt. 
The plots are saved as bandpol_<station>_<exp.code>_<exp.name>.png and 
xcorrmx_<station>_<exp.code>_<exp.name>.png: 
The diagnostic messages are logged in diagnostics_st_IY.txt:

%run select_bandpols.py -s IY -d /data/geodesy/ -o pltIY -x -a
'''



===============================================================================



plot_multcorr.py: 

Generates plots of the band-pols for one experiment on one station,
computing the multiple correlation coefficients of each band-pol with respect
to other band-pols to reject the plots with the mult-corr below the threshold.

Arguments:
  -t <threshold>              for multiple correlation coefficient, 0. to 100.
  -s <a station letter>, like E, G, H ... (or in lower case, e, g, h ...);
  -d <pcc_datfiles directory>       like /data/geodesy/3686/pcc_datfiles
  -o <output directory name>        where .png graphs and .txt logs are saved
  -p                                show plot
  -h                                print this text'''

===============================================================================



plot_multcorr_iter.py:

Generate plots of the band-pols for one experiment on one station.

Arguments:
  -t <threshold>  initial value for multiple correlation coefficient, 0. to 100.
  -s <a station letter>, like E, G, H ... (or in lower case, e, g, h ...);
  -d <pcc_datfiles directory>       like /data/geodesy/3686/pcc_datfiles
  -o <output directory name>        where .png graphs and .txt logs are saved
  -p                                show plot (along with saving in .png file)
  -h                                print this text

Station (-s) and data files directory (-d) are mandatory.
If output directory not specified, the results .png and .txt are saved in 
current directory.

===============================================================================



plot_multcorr_group.py:

Generate plots of all the band-pols for one or all the stations,
computing the multiple correlation coefficients of each band-pol with respect
to other band-pols to reject the plots with the mult-corr below the threshold.

Arguments (optional):
  -s <a station letter>, like E (or in lower case, e);
  -o <output directory name> 

If called with no arguments, it processes all the stations E, G, H, I, V, and Y.

===============================================================================
