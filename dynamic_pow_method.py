# This script fits the dynamic powerlaw method to the burst cooling tails


import os
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
import math
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


# Reading a list of datafiles from prebbfit.txt
# example: 30054-04-02-000_3/analysis/GS1826m24_B00_bbfit_fabs.dat
files=[]
bbfit = open('prebbfit.txt','r')
for line in bbfit :
    files.append(line.rstrip('\n'))
bbfit.close()



# List of doublebursts, relevant only to 4U1636-536
doublebursts = []
pre = open('burst_characteristics.txt','r')
for line in pre :
    ad = line.rstrip('\n').split()
    if len(ad) >= 2:
        if ad[1] == 'double':
            doublebursts.append(ad[0])
    if len(ad) == 3:
        if ad[2] == 'double':
            doublebursts.append(ad[0])
pre.close()


# Creating new output files
fluxoutput = open('flux_xi.txt', 'w')
fluxoutput.write("# Burstid, flux, flux_min, flux_max, alpha, alpha_min, alpha_min\n")
fluxoutput.close()



# Main loop over each datafile (each burst) listed in prebbfit.txt
for datafile in files:

    # Tbb in keV
    tbb=[]
    tbberr_low=[]
    tbberr_high = []
    # K=(Rbb [km]/D_10kpc)^2
    kbb=[]
    kbberr_low=[]
    kbberr_high = []
    kbb4=[]
    kbb4err=[]
    # Calculated flux
    fbb=[]
    fbberr_low=[]
    fbberr_high = []
    # Model flux (cflux)
    fbol=[]
    fbolerr=[]
    time=[]
    timeerr=[]
    
    time0=0.0
    fluxunit = 1e-7
    i=0

    #Reading the data from a file
    f = open(datafile,'r')
    for line2 in f:          
        ab=[]
        ac=[]
        if not line2.startswith("#"):   	# Ignoring comments
            ac=line2.rstrip('\n').split()   # converting line into list, separated by whitespaces
            ab = [float(j) for j in ac]     # converting list of strings into a list of floats
            timeerr.append((ab[0]-time0)/2.0)
            time.append(ab[0])
            time0=ab[0]
            tbb.append(ab[7])
            tbberr_low.append(abs(ab[8]-ab[7]))
            tbberr_high.append(abs(ab[9]-ab[7]))
            kbb.append(ab[10])
            kbberr_low.append(abs(ab[11]-ab[10]))
            kbberr_high.append(abs(ab[12]-ab[10]))
            fbol.append(ab[15]/100.0) # flux in units 10^-7
            fbolerr.append(max(abs(ab[17]-ab[15]),abs(ab[16]-ab[15]))/100.0)
            i=i+1

    f.close()
    timeerr[0]=0.25/2.0
    starttime = time[0]
    timenp = np.array(time)

    # const = sigma * (keV/K)^4 * (10^5)^2/ (3.086*10^22)^2 = 1.0781e-11
    const = 1.0781e-11

    # Calculate flux based on normalisation and temperature
    k=0
    for x in timeerr:
        time[k]=time[k]-starttime  #Scale time to start from zero

        # Flux = (Rbb/D)^2 sigma Tbb^4 = Kbb * tbbkeV^4 * const
        fbb.append(const*kbb[k]*tbb[k]**4/fluxunit)
        fbberr_low.append(const*((kbberr_low[k]*tbb[k]**4)+(4*kbb[k]*tbberr_low[k]*tbb[k]**3))/fluxunit)
        fbberr_high.append(const*((kbberr_high[k]*tbb[k]**4)+(4*kbb[k]*tbberr_high[k]*tbb[k]**3))/fluxunit)
        k=k+1 


    name=datafile.split('/')
    print("\n" + "=" * 20 + "\n" + name[0] + "\n" + "=" * 20)




    # Peak flux
    fbol_peak = max(fbol)
    # fp_i is the index of peak flux
    fp_i = fbol.index(fbol_peak)





    # ======================
    # ======================

    # Fitting and plotting 

    # ======================
    # ======================

    #define pow to be fitted in log scale, i.e. a line
    def func(xdata, a, b):
        return a*xdata + b


    # Initialize plot
    fig = plt.figure()
    plt.rcParams.update({'font.size': 10})
    ax2 = fig.add_subplot(221)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_title('Log-log')
    ax2.minorticks_on()
    ax2.errorbar(timenp, fbb, yerr=[fbberr_low, fbberr_high], fmt='k.')
    ax22 = fig.add_subplot(222, sharex=ax2)
    data = mlines.Line2D([], [], color='k', marker='.', linestyle='solid', label='Flux data with calculated errors')
    fit = mlines.Line2D([], [], color='b', linestyle='solid', label='Fit 0.7-0.05 of peak')
    ax2.legend(handles=[data, fit], prop={'size':6}, loc=1)   


    # Create empty arrays 
    indextimes_pow = []
    xis = []
    chi2Pow = []
    PowPars1 = []
    fluxfit=[]
    fluxerrfit=[]
    timefit=[]
    timefit_pow=[]
    PowErrors_low = []
    PowErrors_high = []

    ###
    ### Fitting one pow to the whole cooling tail
    ###

    # Looping through fluxes and adding to a table the fluxes between 0.7-0.001 of peak flux
    i=0
    for z in fbb:
        if i > fp_i and z < 0.7*fbb[fp_i] and z > 0.001*fbb[fp_i]:
            fluxfit.append(math.log10(z))
            fluxerrfit.append(math.log10(max(fbberr_low[i], fbberr_high[i])))
            if time[i] == 0:
                timefit_pow.append(0.1)
            else:
                timefit_pow.append(math.log10(time[i]))
        i = i + 1

    # Fit a powerlaw to the flux got from previous loop
    x0 = np.array([0.0, 0.0])
    parsP, covarP = opt.curve_fit(func, timefit_pow, fluxfit, x0, fluxerrfit)

    # Creating "array" of the pow index for plotting
    xi1=[]
    for y in timefit:
        xi1.append(-parsP[0])
    # Plotting the results
    timefit = np.array(timefit)
    ax22.plot(timefit, xi1, 'k-', linewidth=1.5)
    ax22.plot(timefit[0], xi1[0], 'k.')
    ax22.plot(timefit[len(timefit)-1], xi1[len(tau1)-1], 'k.')
    ax2.plot(timefit, 10**parsP[1]*timefit**parsP[0], 'b-')

    ###
    ### Dynamic pow fitting    
    ###

    # Looping through the burst flux and fit a pow into the time-window at each time
    for j in range(0, len(fbb)-5):
        
        flux = []
        fluxerr_low = []
        fluxerr_high = []
        fluxfit=[]
        fluxerrfit=[]
        timefit_pow=[]

        # Looping through the flux-array to get the fluxes and times inside the time window
        i=0
        for z in fbb:           
            if i >= j and i < j+7:  # Define window size
                flux.append(z)
                fluxerr_low.append(fbberr_low[i])
                fluxerr_high.append(fbberr_high[i])
                fluxfit.append(math.log10(z))
                fluxerrfit.append(math.log10(max(fbberr_low[i], fbberr_high[i])))
                if time[i] == 0:    #if time is zero, logarithm explodes
                    timefit_pow.append(0.1)
                else:
                    timefit_pow.append(math.log10(time[i]))
            i=i+1

        fluxfit = np.array(fluxfit)
        fluxerrfit = np.array(fluxerrfit)
        timefit_pow = np.array(timefit_pow)

        # Fit pow to flux inside the window
        x0 = np.array([0.0, 0.0])
        parsPow, covarPow = opt.curve_fit(func, timefit_pow, fluxfit, x0, fluxerrfit)

        # Saving results into an array
        indextimes_pow.append(timefit[3])
        xis.append(-parsPow[0])



        #Chi^2 values with asymmetric errors:
        def chi2(par):
            chi2 = 0
            i=0
            for z in flux:
                if z - 10**parsPow[1]*timefit[i]**parsPow[0] < 0:
                    chi2 = chi2 + ((z - 10**parsPow[1]*timefit[i]**par)/fluxerr_high[i])**2
                else:
                    chi2 = chi2 + ((z - 10**parsPow[1]*timefit[i]**par)/fluxerr_low[i])**2
                i = i + 1
            return chi2


        # Parameter errors
        chi2Pow_bestfit = chi2(parsPow[0])
        PowErr_low = 0
        PowErr_high = 0
        for h in range(0, 1000, 1):
            h = h/1000.0
            if PowErr_low == 0 and abs(chi2Pow_bestfit - chi2(parsPow[0] - h)) >= 1:
                PowErr_low = h
            if PowErr_high == 0 and abs(chi2Pow_bestfit - chi2(parsPow[0] + h)) >= 1:
                PowErr_high = h

        # Setting hard limits for the parameter erros, otherwise plots can be very ugly
        if PowErr_low == 0:
            PowErr_low = 1.0
        if PowErr_high == 0:
            PowErr_high = 1.0

        # error values into arrays
        PowErrors_low.append(PowErr_low)
        PowErrors_high.append(PowErr_high)


        # If burst is a double burst, the starting times of the SECOND cooling tails are in a file starttimes2.txt
        # here we check if doubleburst and then read the starting time of the second cooling tail
        # This is relevant only to 4U1636-53
        savetime = 0
        if name[0] in doublebursts:
            stime = open('starttimes2.txt','r')
            for line in stime:
                if not line.startswith("#"):  
                    ad = line.rstrip('\n').split()
                    if ad[0] == name[0]:
                        savetime = float(ad[2])
            stime.close()

        # If double burst, save values into file only from the second cooling tail
        if savetime != 0:
            if timefit[0] >= savetime and parsPow[0] < 0:
                #print "\n! ! ! ! !\n " + name[0] + "\n ! ! ! ! !\n"
                fluxoutput = open('flux_xi.txt', 'a')
                fluxoutput.write(name[0] + " " + str(flux[3]) + " " + str(fluxerr_low[3]) + " " + str(fluxerr_high[3]) + " " + str(parsPow[0]) + " " + str(PowErr_low) + " " + str(PowErr_high) + "\n")
                fluxoutput.close()
        # Otherwise save starting from the peak flux
        else:
            if j >= fp_i:
                fluxoutput = open('flux_xi.txt', 'a')
                fluxoutput.write(name[0] + " " + str(flux[3]) + " " + str(fluxerr_low[3]) + " " + str(fluxerr_high[3]) + " " + str(parsPow[0]) + " " + str(PowErr_low) + " " + str(PowErr_high) + "\n")
                fluxoutput.close()






    #####################################
    # PLOTTING
    #####################################

    xfill = np.linspace(0.1, 100, 1000)
    ax22.set_ylabel('Xi')
    ax22.set_title('Powerlaw index of each 7 point fit')
    ax22.fill_between(xfill, 1.25, 1.33, color='m', alpha=0.2)
    ax22.fill_between(xfill, 1.67, 2.0, color='r', alpha=0.2)
    ax22.errorbar(indextimes_pow, xis, yerr=[PowErrors_low, PowErrors_high], fmt='b.')
    ax22.plot(indextimes_pow, xis, 'b-')
    ax22.minorticks_on()

    for item in ax22.title:
        item.set_fontsize(8)

    ions = mpatches.Patch(color='m', alpha=0.2, label='Ions')
    electrons = mpatches.Patch(color='r', alpha=0.2, label='Electrons')
    wholefit = mlines.Line2D([], [], color='k', marker='.', linestyle='solid', label='Fit 0.7-0.05 of peak')
    xis = mlines.Line2D([], [], color='b', marker='.', linestyle='solid', label='Fit 7 points')
    ax22.legend(handles=[ions, electrons, wholefit, xis], prop={'size':6}, loc=2)    

	#Plot full burst
    ax5 = fig.add_subplot(223)
    ax5.errorbar(timenp, fbol, yerr=[fbberr_low, fbberr_high], fmt='k.')
    ax5.set_title('Whole burst')
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('Flux')
    ax5.minorticks_on()
    ax5.set_xlim(-2, max(time)+2)

	#Plot temperature and normalisation of full burst
    ax7 = fig.add_subplot(326, sharex=ax5)
    ax7.errorbar(timenp, kbb, yerr=[kbberr_low, kbberr_high], color='k')
    ax7.set_ylabel('Normalisation', color='k')
    for tl in ax7.get_yticklabels():
        tl.set_color('k')
    ax7.set_xlim(-2, max(time)+2)
    ax7.minorticks_on()
    ax7.set_xlabel('Time (s)')
    ax6 = ax7.twinx()
    ax6.set_xlim(-2, max(time)+2)
    ax6.errorbar(timenp, tbb, yerr=[tbberr_low, tbberr_high], color='r')
    ax6.set_ylabel('Temperature', color='r')
    for tl in ax6.get_yticklabels():
        tl.set_color('r')
    ax6.minorticks_on()

	#Save figure
    source = os.getcwd().split('/')
    plt.suptitle(source[len(source)-1] + '/' + name[0])    
    fig.set_size_inches(13, 9.0)   
    fig.savefig('plot_linlog/' + name[0] + '.png', dpi=100)   
    plt.close()

