# vgos

Unaided selection of good cable delay data.
The coefficient of multiple correlation of each band-pol with respect 
to other seven is used.
Also, medians for the rows (or columns) of the correlation matrix are calcuated 

===============================================================================



plot_multcorr_all.py:

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
  -t <threshold>              for multiple correlation coefficient, 0. to 1.
  -m <threshold>              for correlation median, -1 to 1.
  -s <a station letter>, like E, G, H ... (or in lower case, e, g, h ...);
  -d <pcc_datfiles directory>       like /data/geodesy/3686/pcc_datfiles
  -o <output directory name>        where .png graphs and .txt logs are saved
  -x                          plot cross-correlation matrix
  -p                          show plot in X-window
  -a                          Make .png plots and .txt files for all available
                              data under directory in -d (like -d /data/geodesy)
                              If one or more stations are given in -s, like
                              -s E or -s VIGH, only data for those stations are
                              plotted and saved in -o directory.
                              If -a is present, -p is ignored (for too many 
                              windows would be open).
  -h                          print this text'''


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