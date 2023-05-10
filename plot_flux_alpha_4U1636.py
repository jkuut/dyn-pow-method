import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import math


# List of pre-bursts
PREbursts=[]
doublebursts=[]
pre = open('burst_characteristics.txt','r')
for line in pre :
    if not line.startswith("#"):
        ad = line.rstrip('\n').split()
        if len(ad) >= 2:
            if ad[1] == 'pre':
                PREbursts.append(ad[0])
            if ad[1] == 'double':
                doublebursts.append(ad[0])
        if len(ad) == 3:
            if ad[2] == 'pre':
                PREbursts.append(ad[0])
            if ad[2] == 'double':
                doublebursts.append(ad[0])
pre.close()


# Initialize plot
fig = plt.figure()
ax1 = fig.add_subplot(221, xscale='log', ylim=(0, 4), xlim=(0.005, 1.0))
ax2 = fig.add_subplot(222, xscale='log', ylim=(0, 4), xlim=(0.005, 1.0))
ax3 = fig.add_subplot(223, xscale='log', ylim=(0, 4), xlim=(0.005, 1.0))
ax4 = fig.add_subplot(224, xscale='log', ylim=(0, 4), xlim=(0.005, 1.0))

ax3.set_xlabel(r'Flux ($10^{-7}$ erg/s/cm$^{2}$)')
ax4.set_xlabel(r'Flux ($10^{-7}$ erg/s/cm$^{2}$)')

ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax4.minorticks_on()

matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 15})


# Set titles to correct places
ax1label = ax1.set_ylabel(r'Powerlaw index $\alpha$')
title = ax1.set_title('Normal bursts')
offset = np.array([-0.2, 0.2])
title.set_position(ax1label.get_position() + offset)
title.set_rotation(90)

ax3label = ax3.set_ylabel(r'Powerlaw index $\alpha$')
title = ax3.set_title('PRE bursts')
offset = np.array([-0.2, 0.1])
title.set_position(ax3label.get_position() + offset)
title.set_rotation(90)

ax2title = ax2.set_title('Hard state')
title = ax4.set_title('Soft state')
offset = np.array([-1.2, 1.2])
title.set_position(ax2title.get_position() + offset)
    




# Find out which burst happened in which state and which is PRE-burst
# burst_hardness.txt should contain a list of burstids and say in which state they happened

hard_norm = []
hard_pre = []
soft_norm = []
soft_pre = []

datafile = 'burst_hardness.txt'
f = open(datafile,'r')
for line2 in f:          
    ac=[]
    if not line2.startswith("#"):
        ac=line2.rstrip('\n').split() 
        if ac[len(ac)-1] == 'hard' and ac[0] in PREbursts and ac[0] not in doublebursts:
            hard_pre.append(ac[0])
        elif ac[len(ac)-1] == 'hard' and ac[0] not in PREbursts and ac[0] not in doublebursts:
            hard_norm.append(ac[0])
        elif ac[len(ac)-1] == 'soft' and ac[0] in PREbursts and ac[0] not in doublebursts:
            soft_pre.append(ac[0])
        elif ac[len(ac)-1] == 'soft' and ac[0] not in PREbursts and ac[0] not in doublebursts:
            soft_norm.append(ac[0])
f.close()


# Read the fluxes and pow indexes of each 4 burst types
flux_hard_norm = []
flux_hard_pre = []
flux_soft_norm = []
flux_soft_pre = []
xi_hard_norm = []
xi_hard_pre = []
xi_soft_norm = []
xi_soft_pre = []
f = open('flux_xi.txt', 'r')
for line in f:
    ad = line.rstrip('\n').split()
    if not line.startswith("#"):
        if -float(ad[7]) < 4.1 and -float(ad[7]) > -0.1: #Ignore outliers for nicer plots
            if ad[0] in hard_norm:
                flux_hard_norm.append(float(ad[1]))
                xi_hard_norm.append(-float(ad[7]))
            elif ad[0] in hard_pre:
                flux_hard_pre.append(float(ad[1]))
                xi_hard_pre.append(-float(ad[7]))
            elif ad[0] in soft_pre:
                flux_soft_pre.append(float(ad[1]))
                xi_soft_pre.append(-float(ad[7]))
            elif ad[0] in soft_norm:
                flux_soft_norm.append(float(ad[1]))
                xi_soft_norm.append(-float(ad[7]))
f.close()

# Lazy way of making sure the bin sizes are the same in each subplot...
# Add new datapoints just outside the plot limits and use those corner points to define the bin sizes
flux_hard_pre.append(1.1)
flux_hard_norm.append(1.1)
flux_soft_norm.append(1.1)
flux_soft_pre.append(1.1)
flux_hard_pre.append(0.0001)
flux_hard_norm.append(0.0001)
flux_soft_norm.append(0.0001)
flux_soft_pre.append(0.0001)
xi_hard_pre.append(4.1)
xi_hard_norm.append(4.1)
xi_soft_norm.append(4.1)
xi_soft_pre.append(4.1)
xi_hard_pre.append(-0.1)
xi_hard_norm.append(-0.1)
xi_soft_norm.append(-0.1)
xi_soft_pre.append(-0.1)

grids=[35,12]


# Plot data
#ax1.hexbin(flux_soft_norm, xi_soft_norm, cmap=plt.cm.gist_heat_r, xscale='log', gridsize=grids)
ax2.hexbin(flux_hard_norm, xi_hard_norm, cmap=plt.cm.gist_heat_r, xscale='log', gridsize=grids)
ax3.hexbin(flux_soft_pre, xi_soft_pre, cmap=plt.cm.gist_heat_r, xscale='log', gridsize=grids)
ax4.hexbin(flux_hard_pre, xi_hard_pre, cmap=plt.cm.gist_heat_r, xscale='log', gridsize=grids)


## Create the colourbar
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
data = ax1.hexbin(flux_soft_norm, xi_soft_norm, cmap=plt.cm.gist_heat_r, xscale='log', gridsize=grids)
cb = plt.colorbar(data, cmap=plt.cm.gist_heat_r, cax=cbar_ax)
cb.set_ticks(np.linspace(data.get_array().min(), data.get_array().max(), 3))
cb.set_ticklabels(np.linspace(0, 1, 3))
cb.set_label('Frequency')



## Define functions that are used in avg and std calculation

# Highest posterior density interval
# Found using iterating the level of some arbitratry distribution where
# we have some fraction of points
def limit_iter(hist, level):
    acc = 0.0001
    diff = 1.0
    iter=0
    sum0 = sum(hist)
    left=0.0
    right=1.0
    midpoint=0.5
    while np.abs(diff) > acc and iter < 200:
        midpoint = (right + left) / 2.0
        upmid = 0.0
        for i in range(len(hist)):
            if hist[i] > midpoint:
                upmid += hist[i]
        ratio = upmid/sum0
        diff = ratio - level
        if diff < 0.0:
            right = midpoint
        else:
            left = midpoint
        iter += 1
    return midpoint
#using previous level finder locate the edges corresponding to these
# limits
def conf_lims(hist, bin_edges, level):
    hist = hist/max(hist)
    bins = bin_edges[0:-1] + np.diff(bin_edges)
    hlim = limit_iter(hist, level)
    i1 = 0
    xlo = bin_edges[i1]
    for j in range(len(bin_edges)):
        if hist[j] >= hlim:
            i1 = j
            xlo = bin_edges[i1]
            break
    i2 = len(bin_edges)-1
    xhi = bin_edges[i2]
    for j in reversed(range(len(bin_edges))):
        if hist[j-1] >= hlim:
            i2 = j
            xhi = bin_edges[i2]
            break
    return xlo, xhi, i1, i2, hlim


# Calculate average and standard deviation 

def avg_and_std(flux, xis, fluxBins = 20, level=0.68):

    levels = np.linspace(min(np.log10(flux)), max(np.log10(flux)), fluxBins)
    avg_xis = []
    std_xis_min = []
    std_xis_max = []
    avg_fluxes = []
    std_xis = []

    for i in range(1, fluxBins):
        j = 0
        xis_temp = []
        for fl in flux:
            if np.log10(fl) >= levels[i-1] and np.log10(fl) < levels[i]:
                xis_temp.append(xis[j])
            j += 1
        if len(xis_temp) > 0:
            avg_xis.append(np.mean(xis_temp))
            avg_fluxes.append(10**((levels[i-1] + levels[i])/2.0))
            hist, bin_edges = np.histogram(xis_temp, bins=50, normed=True)
            xlow, xhigh, i1, i2, hlim = conf_lims(hist, bin_edges, level)
            std_xis_min.append(xlow)
            std_xis_max.append(xhigh)
            std_xis.append(np.std(xis_temp))


    avg_fluxes = np.array(avg_fluxes)
    avg_xis = np.array(avg_xis)
    std_xis_min = np.array(std_xis_min)
    std_xis_max = np.array(std_xis_max)

    return avg_fluxes, avg_xis, std_xis_min, std_xis_max, std_xis


lwd = 3.0
alph = 0.7

# Calculate and then plot the avg and std in each subplot
avg_flux, avg_xi, std_xi_min, std_xi_max, std_xis = avg_and_std(flux_soft_norm, xi_soft_norm)
ax1.plot(avg_flux, avg_xi, 'b-', linewidth=lwd)
ax1.plot(avg_flux, std_xi_min, 'b--', alpha=alph, linewidth=lwd)
ax1.plot(avg_flux, std_xi_max, 'b--', alpha=alph, linewidth=lwd)

avg_flux, avg_xi, std_xi_min, std_xi_max, std_xis = avg_and_std(flux_hard_norm, xi_hard_norm, 25)
ax2.plot(avg_flux, avg_xi, 'b-', linewidth=lwd)
ax2.plot(avg_flux, std_xi_min, 'b--', alpha=alph, linewidth=lwd)
ax2.plot(avg_flux, std_xi_max, 'b--', alpha=alph, linewidth=lwd)

avg_flux, avg_xi, std_xi_min, std_xi_max, std_xis = avg_and_std(flux_soft_pre, xi_soft_pre, 25)
ax3.plot(avg_flux, avg_xi, 'b-', linewidth=lwd)
ax3.plot(avg_flux, std_xi_min, 'b--', alpha=alph, linewidth=lwd)
ax3.plot(avg_flux, std_xi_max, 'b--', alpha=alph, linewidth=lwd)

avg_flux, avg_xi, std_xi_min, std_xi_max, std_xis = avg_and_std(flux_hard_pre, xi_hard_pre, 15)
ax4.plot(avg_flux, avg_xi, 'b-', linewidth=lwd)
ax4.plot(avg_flux, std_xi_min, 'b--', alpha=alph, linewidth=lwd)
ax4.plot(avg_flux, std_xi_max, 'b--', alpha=alph, linewidth=lwd)



#Function to plot the cooling models of Cumming et al. (2004) and in't Zand et al. (2014)
def simple_models():
    flux_10_8 = []
    flux_10_9 = []
    flux_10_10 = []
    xi_10_8 = []
    xi_10_9 = []
    xi_10_10 = []
    flux_10_8_onezone = []
    xi_10_8_onezone = []

    fnorm=5.0*10**(-26)  # Normalise the flux, this value is different for each source

    # Read the data of the models
    f1 = open('10_8.csv', 'r')
    for line in f1:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_8.append(fnorm*float(ad[0]))
            xi_10_8.append(float(ad[1]))
    f1.close()
    f2 = open('10_9.csv', 'r')
    for line in f2:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_9.append(fnorm*float(ad[0]))
            xi_10_9.append(float(ad[1]))
    f2.close()
    f3 = open('10_10.csv', 'r')
    for line in f3:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_10.append(fnorm*float(ad[0]))
            xi_10_10.append(float(ad[1]))
    f3.close()
    f4 = open('10_8_onezone.csv', 'r')
    for line in f4:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_8_onezone.append(fnorm*float(ad[0]))
            xi_10_8_onezone.append(float(ad[1]))
    f4.close()

    lwd=2.0
    aa,=ax1.plot(flux_10_8, xi_10_8, 'c-', linewidth=lwd)
    bb,=ax1.plot(flux_10_9, xi_10_9, 'm-', linewidth=lwd)
    cc,=ax1.plot(flux_10_10, xi_10_10, 'y-', linewidth=lwd)
    ax2.plot(flux_10_8, xi_10_8, 'c-', linewidth=lwd)
    ax2.plot(flux_10_9, xi_10_9, 'm-', linewidth=lwd)
    ax2.plot(flux_10_10, xi_10_10, 'y-', linewidth=lwd)
    ax3.plot(flux_10_8, xi_10_8, 'c-', linewidth=lwd)
    ax3.plot(flux_10_9, xi_10_9, 'm-', linewidth=lwd)
    ax3.plot(flux_10_10, xi_10_10, 'y-', linewidth=lwd)
    ax4.plot(flux_10_8, xi_10_8, 'c-', linewidth=lwd)
    ax4.plot(flux_10_9, xi_10_9, 'm-', linewidth=lwd)
    ax4.plot(flux_10_10, xi_10_10, 'y-', linewidth=lwd)
    onezone,=ax1.plot(flux_10_8_onezone, xi_10_8_onezone, 'g-', linewidth=lwd)
    ax2.plot(flux_10_8_onezone, xi_10_8_onezone, 'g-', linewidth=lwd)
    ax3.plot(flux_10_8_onezone, xi_10_8_onezone, 'g-', linewidth=lwd)
    ax4.plot(flux_10_8_onezone, xi_10_8_onezone, 'g-', linewidth=lwd)
    ax2.legend([aa, bb, cc, onezone], [r'$10^{8}$ g/cm$^2$', r'$10^{9}$ g/cm$^2$', r'$10^{10}$ g/cm$^2$', r'Onezone, $10^{8}$ g/cm$^2$'], prop={'size':10})  


# Function to plot some example Kepler models.
# Not really needed/used anymore
def kepler():

    lwd=2.0
    fnorm=1.3*10**(-39)
    j=0
    f1 = open('/scratch/database/db/Kepler/kepler_bursts.txt', 'r')
    for line in f1:
        flux_model = []
        xi_model = []
        if not line.startswith("#"):
            ad = line.rstrip('\n').split()
            f2 = open(ad[len(ad)-1], 'r')
            for line2 in f2:
                if not line2.startswith("#"):
                    ac = line2.rstrip('\n').split()
                    if ad[0] == ac[0]:
                        flux_model.append(fnorm*float(ac[2]))
                        xi_model.append(-float(ac[5]))
            if ad[1] == '20.8':
                mrk = 'o'
                col = 'g'
            elif ad[1] == '48.4':
                mrk = 's'
                col = 'c'
            elif ad[1] == '64.96':
                mrk = '^'
                col = 'b'
            elif ad[1] == '70.48':
                mrk = 'D'
                col = 'k'
            ax1.plot(flux_model, xi_model, color=col, linewidth=lwd, alpha=0.7)
            ax2.plot(flux_model, xi_model, color=col, linewidth=lwd, alpha=0.7)
            ax3.plot(flux_model, xi_model, color=col, linewidth=lwd, alpha=0.7)
            ax4.plot(flux_model, xi_model, color=col, linewidth=lwd, alpha=0.7)
            j+=1
            f2.close()
    f1.close()
    H208 = mlines.Line2D([], [], color='g', label='X=0.208')
    H484 = mlines.Line2D([], [], color='c', label='X=0.484')
    H6496 = mlines.Line2D([], [], color='b', label='X=0.6496')
    H7048 = mlines.Line2D([], [], color='k', label='X=0.7048')
    ax1.legend(handles=[H208, H484, H6496, H7048], prop={'size':8}, loc=2)   


# Plot kepler or simple models?
#kepler()
simple_models()

#Save figure
fig.set_size_inches(13, 9)
fig.savefig('1636_flux_alpha.pdf', dpi=200)
plt.show()

    
    







