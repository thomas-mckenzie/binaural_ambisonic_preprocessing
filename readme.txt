This folder should include all the files necessary to run the Ambisonic HRTF pre-processing code. 

See \test_scripts folder for four test scripts. 

Run load_ambisonic_configuration to generate binaural Ambisonic decoders with the pre-processing techniques. Configure this script to select the desired pre-processing techniques. It will produce one pre-processed binaural Ambisonic decoder, and one with no pre-processing (labelled NPP). The included HRTFs are from the Benjamin Bernschutz KU 100 database: http://audiogroup.web.th-koeln.de/ku100hrir.html

To test the decoders, run test_ambisonic_decoder. This will produce a figure with the diffuse-field responses of both decoders. It will also compare binaural Ambisonic rendering using the two decoders with standard (non-Ambisonic) HRIRs in the folder \hrirs. This comparison is done for perceptual spectral difference (PSD), interaural level difference (ILD) and interaural time difference (ITD). These will be plotted too. 

To listen to the differences between the decoders, run test_sound_percussion, to hear the standard binaural Ambisonic decoder, followed by the pre-processed one, and then an HRTF convolution render. The sound is several percussive sounds virtually panned to various locations (see thesis for exact locations). 
For another test sound script to hear the differences between the decoders, run test_sound_soundscape to hear the standard binaural Ambisonic decoder followed by the pre-processed one. The soundscape excerpts are from the Marc Green EigenScape database: https://zenodo.org/record/1012809

Make sure to de-zip the hrir folder!

Thomas McKenzie, University of York, 2019. Now at University of Edinburgh.
thomas.mckenzie@ed.ac.uk

----------
This code requires the Matlab toolboxes Phased Array System Toolbox and Curve Fitting Toolbox to run. 

This work uses functions written by others. In all cases, these are documented in the code when used. Some functions have been adapted.
This code utilises Ambisonic functions by Archontis Politis: https://uk.mathworks.com/matlabcentral/fileexchange/54833-higher-order-ambisonics-hoa-library
inverse filtering code by Matthes: https://uk.mathworks.com/matlabcentral/fileexchange/19294-inverse-fir-filter
voronoi sphere code from Bruno Luong: https://uk.mathworks.com/matlabcentral/fileexchange/40989-voronoi-sphere
quadrature code from John Burkardt: https://people.sc.fsu.edu/~jburkardt/datasets/sphere_grid/sphere_grid.html
iso226 code by Jeff Tackett: https://uk.mathworks.com/matlabcentral/fileexchange/7028-iso-226-equal-loudness-level-contour-signal?focused=22bd2900-9b20-b731-8ea3-64189da018b4&tab=function
ERB code from IOSR toolbox (Surrey University): https://github.com/IoSR-Surrey/MatlabToolbox
