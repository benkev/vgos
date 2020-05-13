# vgos

Unaided selection of good cable delay data.
The coefficient of multiple correlation of each band-pol with respect 
to other seven is used.



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