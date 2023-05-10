Scripts for the flux-alpha project

The burst data should be in the same format as produced by xspec_burst_fit_v7.pl

dynamic_pow_method.py 
    This script reads the data for each burst listed in prebbfit.txt and then fits a powerlaw to the flux in moving time window.
    The fit results are saved into file flux_xi.txt, and each burst is also plotted separately. 
    
    
plot_flux_alpha_4U1636.py (or plot_flux_alpha_GS1826.py)
    This script reads the fit results from flux_xi.txt and then plots it in a nice way. For sources like 4U1636-53 the bursts
    are separated based on the state and whether PRE happened or not, but for sources like GS1826-24 just one plot is created.
    (That is why two separate scripts)
    Also the simple cooling models of Cumming et al (2004) and in't Zand et al (2014) can be plotted with the data.
    
    
colcol_rms.py
    This script creates the color-color and rms-countrate plots used for state separation in the burst cooling paper. 
