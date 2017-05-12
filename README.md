# CPD
Curie Point Depth Functions.

A set of functions to split up a large aeromag grid (splitgrid) into smaller windows for calculating fft's (runfft),
which are used to calculate Curie Point Depths (calccpd) at each window center.  See 

A. Tanaka, Y. Okubo, O. Matsubayashi
Curie point depth based on spectrum analysis of the magnetic anomaly data in East and Southeast Asia
Tectonophysics, 306 (1999), pp. 461–470

Requires GMT installed for calling grdfft
http://gmt.soest.hawaii.edu/
