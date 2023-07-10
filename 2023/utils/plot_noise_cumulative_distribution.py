import os
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


# Parameters
fullsky = 4 * np.pi * (360/2/np.pi) ** 2 # deg2
nside = 1024
# uK.pixel-1 to uK.arcmin
fac = np.sqrt((360 * 60) ** 2 / np.pi / hp.nside2npix(nside))


# Read Noise variance map
freqs = ['95', '150']
fnames = [ os.path.join(os.environ['ALILENS'], 'Noise_ALI_IHEP_20200730/I_NOISE_%s_C_1024.fits' % freq) for freq in freqs]
# NOTE fits skymap in RING ordering
nlevs = [ hp.read_map(fname, verbose=False, nest=False) * fac for fname in fnames]
nlevs = [ np.array([ 0. if np.isnan(p) else p for p in nlev ])  for nlev in nlevs ]


# Print AliCPT noise map basic information
for _,freq in enumerate(freqs):
    nlev = nlevs[_]
    fsky = len(nlev[nlev>0.]) / len(nlev)
    nlev = nlev[nlev>0.]
    print("Channel %s GHz Noise Map Info" % freq)
    print("    fsky         = %.2f%%" % (fsky*100) )
    print("    nlev min     = %.2f uK.arcmin" % nlev.min())
    print("    nlev max     = %.2f uK.arcmin" % nlev.max())
    print("    nlev mean    = %.2f uK.arcmin" % nlev.mean())
    print("    nlev median  = %.2f uK.arcmin" % np.median(nlev))



# Plot noise cumulative distribution
fig, ax = plt.subplots()

for i,freq in enumerate(freqs):
    nlev = nlevs[i]
    thresholds = np.logspace(np.log10(nlev[nlev!=0.].min()+0.01), np.log10(nlev.max()), 100)

    # Calc cumulative area with unit deg2
    skyarea = lambda threshold : len(nlev[ np.where((nlev>0.)&(nlev<=threshold)) ]) / len(nlev) * fullsky
    # Calc average noise in this area
    nlev_avg = lambda threshold : np.mean( nlev[ np.where((nlev>0.)&(nlev<=threshold) ) ] )

    xs = np.zeros_like(thresholds)
    ys = np.zeros_like(thresholds)
    for _,threshold in enumerate(thresholds):
        xs[_] = skyarea(threshold)
        ys[_] = nlev_avg(threshold)

    ax.plot( xs, ys, label='%s GHz' % freq )


ax.legend(loc='upper left')
ax.set_xlabel('Cumulative Area')
ax.set_ylabel('White Noise RMS [$\mu$K-acrmin]')
ax.set_xlim([0, fullsky * 0.14])
ax.set_ylim([0, 50])
ax.set_title('Map noise cumulative distribution')
ax.grid(which='both', linestyle=':')
fig.savefig('noise_cumulative_distribution.png', dpi=300)
