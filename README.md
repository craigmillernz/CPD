# CPD
Curie Point Depth Functions.

A set of functions to split up a large aeromag grid (splitgrid) into smaller windows for calculating fft's (runfft),
which are used to calculate Curie Point Depths (calccpd) at each window center.  See 

A. Tanaka, Y. Okubo, O. Matsubayashi
Curie point depth based on spectrum analysis of the magnetic anomaly data in East and Southeast Asia
Tectonophysics, 306 (1999), pp. 461–470

A. R. Bansal, G. Gabriel, V. P. Dimri, and C. M. Krawczyk, (2011), "Estimation of depth to the bottom of magnetic sources by a modified centroid method for fractal distribution of sources: An application to aeromagnetic data in Germany," GEOPHYSICS 76: L11-L22.
https://doi.org/10.1190/1.3560017

Change to using pygmt for calling gmt functions directly within python, not as external system calls.

Requires pygmt 
https://www.pygmt.org/dev/overview.html
