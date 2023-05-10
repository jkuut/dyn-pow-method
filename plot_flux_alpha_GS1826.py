import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches



# Initialize plot
fig = plt.figure()
ax1 = fig.add_subplot(111, xscale='log', ylim=(0, 5), xlim=(0.003, 0.4))
plt.rcParams.update({'font.size': 14})
ax1.set_xlabel(r'Flux ($10^{-7}$ erg/s/cm$^{2}$)')
ax1.set_ylabel(r'Powerlaw index $\alpha$')
ax1.minorticks_on()
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 15})



# Read the flux and alpha values
flux = []
xi = []
f = open('flux_xi.txt', 'r')
for line in f:
    ad = line.rstrip('\n').split()
    if not line.startswith("#"):
        if -float(ad[4]) < 5.1 and -float(ad[4]) > -0.1:
            flux.append(float(ad[1]))
            xi.append(-float(ad[4]))
f.close()


#flux.append(1.1)
#xi.append(5.1)
#flux.append(0.0001)
#xi.append(-0.1)
# gridsize=[35,12]

# Plot data and colorbar
data = ax1.hexbin(flux, xi, cmap=plt.cm.gist_heat_r, xscale='log')
cb = plt.colorbar(data, cmap=plt.cm.gist_heat_r)
cb.set_ticks(np.linspace(data.get_array().min(), data.get_array().max(), 3))
cb.set_ticklabels(np.linspace(0, 1, 3))
cb.set_label('Frequency')


# Functions used to calculate the avg and std
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
#using previous level finder locate the edges corresponding to these limits
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



# Calculate avg and std
def avg_and_std(flux, xis, fluxBins = 25, level=0.68):

    levels = np.linspace(min(np.log10(flux)), max(np.log10(flux)), fluxBins)
    avg_xis = []
    std_xis_min = []
    std_xis_max = []
    avg_fluxes = []
    std_xis = []

    for i in range(0, fluxBins):
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
alph = 0.9
# Plot avg and std
avg_flux, avg_xi, std_xi_min, std_xi_max, std_xis = avg_and_std(flux, xi)
ax1.plot(avg_flux, avg_xi, 'b-', linewidth=lwd)
ax1.plot(avg_flux, std_xi_min, 'b--', alpha=0.7, linewidth=lwd)
ax1.plot(avg_flux, std_xi_max, 'b--', alpha=0.7, linewidth=lwd)

# Save the calculated average to a file for later use
f = open('avg_data.txt', 'w')
for i in range(0, len(avg_flux)):
    f.write(str(avg_flux[i]) + ' ' + str(avg_xi[i]) + ' ' + str(std_xi_min[i]) + ' ' + str(std_xi_max[i]) + '\n')
f.close() 


#Plot simple cooling models
def simple_models():
    flux_10_8 = []
    flux_10_9 = []
    flux_10_10 = []
    xi_10_8 = []
    xi_10_9 = []
    xi_10_10 = []
    flux_10_8_onezone = []
    xi_10_8_onezone = []

    fnorm=3.0*10**(-26)
    f1 = open('/scratch/database/db/10_8.csv', 'r')
    for line in f1:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_8.append(fnorm*float(ad[0]))
            xi_10_8.append(float(ad[1]))
    f1.close()
    f2 = open('/scratch/database/db/10_9.csv', 'r')
    for line in f2:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_9.append(fnorm*float(ad[0]))
            xi_10_9.append(float(ad[1]))
    f2.close()
    f3 = open('/scratch/database/db/10_10.csv', 'r')
    for line in f3:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_10.append(fnorm*float(ad[0]))
            xi_10_10.append(float(ad[1]))
    f3.close()
    f4 = open('/scratch/database/db/10_8_onezone.csv', 'r')
    for line in f4:
        ad = line.rstrip('\n').split(',')
        if not line.startswith("#"):
            flux_10_8_onezone.append(fnorm*float(ad[0]))
            xi_10_8_onezone.append(float(ad[1]))
    f4.close()

    lwd=2.5
    aa,=ax1.plot(flux_10_8, xi_10_8, 'c-', linewidth=lwd)
    bb,=ax1.plot(flux_10_9, xi_10_9, 'm-', linewidth=lwd)
    cc,=ax1.plot(flux_10_10, xi_10_10, 'y-', linewidth=lwd)
    onezone,=ax1.plot(flux_10_8_onezone, xi_10_8_onezone, 'g-', linewidth=lwd)
    ax1.legend([aa, bb, cc, onezone], [r'$10^{8}$ g/cm$^2$', r'$10^{9}$ g/cm$^2$', r'$10^{10}$ g/cm$^2$', r'Onezone, $10^{8}$ g/cm$^2$'], prop={'size':18})  

#Plot example kepler models
def kepler():

    lwd=2.5
    fnorm=1.5*10**(-39)
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
            j+=1
            f2.close()
    f1.close()


    H208 = mlines.Line2D([], [], color='g', label='X=0.208')
    H484 = mlines.Line2D([], [], color='c', label='X=0.484')
    H6496 = mlines.Line2D([], [], color='b', label='X=0.6496')
    H7048 = mlines.Line2D([], [], color='k', label='X=0.7048')
    ax1.legend(handles=[H208, H484, H6496, H7048], prop={'size':12}, loc=2)   


#kepler()
simple_models()

plt.rcParams.update({'font.size': 24})
fig.set_size_inches(13, 9)
fig.savefig('pdfplots/1826_combined_cb.pdf', dpi=200)
plt.show()










