# This script creates the color-color and rms plots used in state separation
# This is quite messy because of the different ways the rms and coco files are defined 


import os
import numpy as np
from matplotlib import pyplot as plt
from math import exp
from math import sqrt
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.lines as mlines



# List of pre-bursts
PREbursts=[]
doublebursts=[]
PREbursts2=[]
pre = open('burst_characteristics.txt','r')
for line in pre :
    if not line.startswith("#"):
        ad = line.rstrip('\n').split()
        ax = ad[0].split('_')
        if len(ad) >= 2:
            if ad[1] == 'pre':
                PREbursts.append(ad[0])
                PREbursts2.append(ax[0])
            if ad[1] == 'double':
                doublebursts.append(ad[0])
        if len(ad) == 3:
            if ad[2] == 'pre':
                PREbursts.append(ad[0])
                PREbursts2.append(ax[0])
            if ad[2] == 'double':
                doublebursts.append(ad[0])
pre.close()


# Read the touchdown fluxes of prebursts and calculate the average of those
# This is used as an eddington flux
tdfluxes = []
f = open('burst_properties.txt','r')
for line in f:
    if not line.startswith("#"):
        ae = line.rstrip('\n').split()
        if ae[0] in PREbursts:
            tdfluxes.append(float(ae[7]))
f.close()
tdflux = 0
for j in tdfluxes:
    tdflux = tdflux + j
tdflux = 10**(-7)*tdflux/len(tdfluxes)


# Create lists of hard and soft bursts
# hard and soft lists contain the burstids in form 10088-01-08-01_3
# while hard2 and soft2 contain the burstids in form 10088-01-08-01
hard = []
soft = []
hard2 = []
soft2 = []
f = open('burst_hardness.txt','r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split() 
        if ac[len(ac)-1] == 'hard':
            hard.append(ac[0])
            ax=ac[0].split('_')
            hard2.append(ax[0])
        elif ac[len(ac)-1] == 'soft':
            soft.append(ac[0])
            ax=ac[0].split('_')
            soft2.append(ax[0])
f.close()


# Read the fluxes in four different energy ranges
flux1 = []
flux2 = []
flux3 = []
flux4 = []
flux_total = []
f = open('1636_coco_nobursts.dat','r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split()
        flux1.append(float(ac[1]))
        flux2.append(float(ac[3]))
        flux3.append(float(ac[5]))
        flux4.append(float(ac[7]))
        flux_total.append(float(ac[9]))
f.close()
flux1 = np.array([float(j) for j in flux1])
flux2 = np.array([float(j) for j in flux2])
flux3 = np.array([float(j) for j in flux3])
flux4 = np.array([float(j) for j in flux4])

#Persistent flux divided by eddington flux
flux_total = np.array([float(j) for j in flux_total])/tdflux
# Hard and soft colours
hard_all = flux4/flux3
soft_all = flux2/flux1


# Find out the hard and soft colours and the persistent flux of each burst
hard_burst_hard = []
hard_burst_soft = []
soft_burst_hard = []
soft_burst_soft = []
hard_burst_hard_pre = []
hard_burst_soft_pre = []
soft_burst_hard_pre = []
soft_burst_soft_pre = []
fluxper_hard = []
fluxper_soft = []
fluxper_hard_pre = []
fluxper_soft_pre = []
f = open('burst_hardness.txt','r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split()
        if ac[0] not in PREbursts:
            if ac[0] in hard:
                hard_burst_hard.append(float(ac[4]))
                hard_burst_soft.append(float(ac[6]))
                fluxper_hard.append(float(ac[2])/tdflux)
            elif ac[0] in soft:
                soft_burst_hard.append(float(ac[4]))
                soft_burst_soft.append(float(ac[6]))
                fluxper_soft.append(float(ac[2])/tdflux)
        elif ac[0] in PREbursts:
            if ac[0] in hard:
                hard_burst_hard_pre.append(float(ac[4]))
                hard_burst_soft_pre.append(float(ac[6]))
                fluxper_hard_pre.append(float(ac[2])/tdflux)
            elif ac[0] in soft:
                soft_burst_hard_pre.append(float(ac[4]))
                soft_burst_soft_pre.append(float(ac[6]))
                fluxper_soft_pre.append(float(ac[2])/tdflux)
f.close()


# Plotting
fig = plt.figure()
ax = fig.add_subplot(131)
ax.minorticks_on()
ax.plot(soft_all, hard_all, color='0.75',marker='.',linestyle='none')
#ax.set_xlabel('Soft color (4$-$6.4 keV)/(3$-$4 keV)')
ax.set_ylabel('Hard color (9.7$-$16 keV)/(6.4$-$9.7 keV)')

rmshard,=ax.plot(hard_burst_soft, hard_burst_hard, 'ko')
rmssoft,=ax.plot(soft_burst_soft, soft_burst_hard, 'bo')
rmshard_pre,=ax.plot(hard_burst_soft_pre, hard_burst_hard_pre, 'k^', markeredgecolor='grey')
rmssoft_pre,=ax.plot(soft_burst_soft_pre, soft_burst_hard_pre, 'b^', markeredgecolor='c')

ax3 = fig.add_subplot(132)
ax3.minorticks_on()
ax3.set_xscale('log')
ax3.set_xlim(0.01, 0.5)
ax3.plot(flux_total, hard_all, color='0.75',marker='.',linestyle='none')
#ax3.set_xlabel(r'Persistent flux F$_{\mathrm{per}}$/<F$_{\mathrm{td}}$>')
ax3.set_ylabel('Hard color (9.7$-$16 keV)/(6.4$-$9.7 keV)')
ax3.plot(fluxper_hard, hard_burst_hard, 'ko')
ax3.plot(fluxper_soft, soft_burst_hard, 'bo')
ax3.plot(fluxper_hard_pre, hard_burst_hard_pre, 'k^', markeredgecolor='grey')
ax3.plot(fluxper_soft_pre, soft_burst_hard_pre, 'b^', markeredgecolor='c')




# Read the rms-values for each burst
rms_hard = []
rms_soft = []
rms_hard_error = []
rms_soft_error = []
rms_hard_pre = []
rms_soft_pre = []
rms_hard_error_pre = []
rms_soft_error_pre = []
f = open('4U1636_rms.dat','r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split() 
        if ac[1] not in PREbursts2:
            if ac[1] in hard2:
                rms_hard.append(float(ac[2]))
                rms_hard_error.append(float(ac[3]))
            elif ac[1] in soft2:
                rms_soft.append(float(ac[2]))
                rms_soft_error.append(float(ac[3]))
        elif ac[1] in PREbursts2:
            if ac[1] in hard2:
                rms_hard_pre.append(float(ac[2]))
                rms_hard_error_pre.append(float(ac[3]))
            elif ac[1] in soft2:
                rms_soft_pre.append(float(ac[2]))
                rms_soft_error_pre.append(float(ac[3]))
f.close()

# Read the count rate (not hardness) for each burst
hardness_hard = []
hardness_soft = []
hardness_hard_error = []
hardness_soft_error = []
hardness_hard_pre = []
hardness_soft_pre = []
hardness_hard_error_pre = []
hardness_soft_error_pre = []
f = open('4U1636_rms.dat','r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split() 
        if ac[1] not in PREbursts2:
            if ac[1] in hard2:
                hardness_hard.append(float(ac[6]))
                hardness_hard_error.append(float(ac[7]))
            elif ac[1] in soft2:
                hardness_soft.append(float(ac[6]))
                hardness_soft_error.append(float(ac[7]))
        elif ac[1] in PREbursts2:
            if ac[1] in hard2:
                hardness_hard_pre.append(float(ac[6]))
                hardness_hard_error_pre.append(float(ac[7]))
            elif ac[1] in soft2:
                hardness_soft_pre.append(float(ac[6]))
                hardness_soft_error_pre.append(float(ac[7]))
f.close()

# Plot count rate vs rms 
ax2 = fig.add_subplot(133)
ax2.minorticks_on()
ax2.plot(hardness_hard, rms_hard, marker='o', color='k', linestyle='none')
ax2.plot(hardness_soft, rms_soft, marker='o', color='b', linestyle='none')
ax2.plot(hardness_hard_pre, rms_hard_pre, marker='^', color='k', markeredgecolor='grey', linestyle='none')
ax2.plot(hardness_soft_pre, rms_soft_pre, marker='^', color='b', markeredgecolor='c', linestyle='none')
ax2.set_ylabel(r'Fractional rms')
ax2.set_ylim(0, 0.23)

# Create legend
hards = mpatches.Patch(color='k', label='Hard bursts')
softs = mpatches.Patch(color='b', label='Soft bursts')
pre = mlines.Line2D([], [], color='k', marker='^', linestyle='none', label='PRE bursts')
norm = mlines.Line2D([], [], color='k', marker='o', linestyle='none', label='Normal bursts')
ax2.legend(handles=[hards, softs, pre, norm], prop={'size':9}, loc=1)

# Typical error bars
ax.errorbar(2.0, 0.65, xerr=0.03, yerr=0.008, fmt='.', color='g')
ax3.errorbar(0.3, 0.65, xerr=0.03, yerr=0.008, fmt='.', color='g')
ax2.errorbar(50, 0.02, xerr=15, yerr=0.005, fmt='.', color='g')

# Set every other tick label invisible 
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
for label in ax2.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)



plt.suptitle('4U 1636-536')

matplotlib.rcParams['pdf.fonttype'] = 42
plt.subplots_adjust(wspace=0.3)
fig.set_size_inches(13.0, 4.0)   
fig.savefig('pdfplots/1636_colcol_rms.pdf', bbox_inches='tight', dpi=200)   
plt.close()
