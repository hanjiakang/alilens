""" ALL Plot Scripts

    * recon_cl plot
    * SNR plot
    * reconstruction map

"""


import os
import sys
import numpy as np
import healpy as hp

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle # Error boxes
from matplotlib.collections import PatchCollection # Error boxes
from plancklens import utils
from plancklens.qresp import get_response

sys.path.insert(0, './')
from one import *
import params as par
import bandpowers


# Error boxes
def make_error_boxes(ax, xdata, ydata, xerr, yerr, facecolor, edgecolor='None', alpha=0.5):
    errorboxes = []
    for x, y, xe, ye in zip(xdata, ydata, xerr.T, yerr.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)





# Parameters
#savePath = './'
#qe_keys = {
#           'p_eb':['k', 'EB'],
#           'p_te':['r', 'TE']
#           } # MV estimator
#btype = 'agr2'

# Parameters
savePath = './'
qe_keys = {
            'ptt':['royalblue', 'TT'] # Temperature only
           ,'p_p':['green', 'Pol only'] # Polarization only
           ,'p':['tomato', 'MV'],
#          'pee':['royalblue', 'EE']
 #          ,'p_p':['green', 'Pol']
 #          ,'p_eb':['tomato', 'EB']
          } # MV estimator
btype = 'agr2'



# Plot Reconstructed power spectrum
fig, ax = plt.subplots()

for key in qe_keys.keys():
    print(key)
    bp = bandpowers.ffp10_binner(key, key, par, btype, ksource='p')
    bells = bp.bin_lavs
    bpower = (bp.get_dat_bandpowers() - bp.get_rdn0() - bp.get_n1()) * bp.get_bmmc()
#    bpower = (bp.get_dat_bandpowers() - bp.get_rdn0())
    
    # Potential Estimator
    yerr = np.sqrt(bp.get_cov().diagonal())
    l, u, c = bandpowers.get_blbubc(btype)
    xerr = np.stack([bp.bin_lavs - l, u - bp.bin_lavs + 1])
    make_error_boxes(ax, bells, bpower, xerr, np.stack([yerr, yerr]), facecolor=qe_keys[key][0], edgecolor='None', alpha=0.2)
    ax.scatter(bells, bpower, label='$C_L^{\hat\phi\hat\phi}-\hat N_L^{(0)}$  (%s)' % qe_keys[key][1], c=qe_keys[key][0], s=13)

# fiducial lensing power spectrum
ax.plot(np.arange(2049), bp.clkk_fid, label=r'$C_L^{\phi\phi, fid}$', c='k')
ax.semilogx()
ax.semilogy()
ax.set_title('Lensing Reconstruction', fontsize=12)
ax.set_xlabel(r'$L$', fontsize=9)
ax.set_ylabel(r'$10^7 L^2 (L + 1)^2 C_L^{\phi\phi}/2\pi$', fontsize=9)
ax.legend(loc='upper right', fontsize=7)
ax.set_xlim(l.min(), u.max())
ax.set_ylim(0.1, 3.)
fig.savefig(os.path.join(ALILENS, 'recon_cl.png'), dpi=300)
#fig.savefig('recon_cl.png', dpi=300)





# Plot SNR plot
fig, ax = plt.subplots()
for qe_key in qe_keys.keys():
    print(qe_key)
    bp = bandpowers.ffp10_binner(qe_key, qe_key, par, btype, ksource='p')
    cov = bp.get_cov()
    signal = (bp.get_dat_bandpowers() - bp.get_rdn0() - bp.get_n1()) * bp.get_bmmc()
#    signal = bp.get_fid_bandpowers()
    l, u, c = bandpowers.get_blbubc(btype)
    SNRs = signal / np.sqrt(cov.diagonal())

    ax.plot(bp.bin_lavs, SNRs, c=qe_keys[qe_key][0], label='SNR (%s)' % qe_key)
    print(np.dot(np.dot(signal, np.matrix(cov).I), signal)[0,0] ** 0.5)
    print(np.sqrt(sum(signal**2 / cov.diagonal())))
#    print(cov)
ax.semilogx()
ax.set_title('Reconstruction SNR', fontsize=12)
ax.set_xlabel(r'$L$', fontsize=12)
ax.set_ylabel('SNR', fontsize=12)
ax.legend()
ax.set_xlim(l.min(), 2048)

fig.savefig(os.path.join(ALILENS, 'recon_snr.png'), dpi=300)







# Reconstruction map
# def view_map(m, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
#     """ View map.
#     """
#     # TODO beautify this plot
#     rot = [180, 60, 0]


#     m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
#     m[ m==0. ] = np.nan # in case the input map is an apodization mask

#     if min==None: min = m[ ~np.isnan(m) ].min()
#     if max==None: max = m[ ~np.isnan(m) ].max()

#     hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.graticule()
#     plt.savefig(savePath, dpi=300)


# mask_b = np.where(hp.read_map(os.path.join(ALILENS, 'sims/ninv/ninv_t.fits')) > 0, 1, 0)
# lmax = 2048
# q2k = lambda l: l*(l + 1) / 2 # potential -> convergence
# q2d = lambda l: (l*(l + 1)) ** 0.5 # potential -> deflection
# cut = np.where((np.arange(lmax + 1) > 8) * (np.arange(lmax + 1) < 2000), 1, 0) # band limit

# # wiener filter
# wiener_dat = np.loadtxt(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/nlkk.dat')).transpose()
# wiener = (wiener_dat[2] - wiener_dat[1]) * utils.cli(wiener_dat[2])

# # input deflection map
# qlm_input = hp.map2alm(hp.read_map(os.path.join(ALILENS, 'sims/cmbs/map_P_1024_0300.fits')))
# dlm_input = hp.almxfl(qlm_input, cut * q2d(np.arange(lmax + 1)))
# dmap_input = hp.alm2map(dlm_input, nside=1024)

# # reconstruction map
# klm_recon = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/dat_klm.fits'))
# dlm_recon = hp.almxfl(klm_recon, cut * utils.cli(q2k(np.arange(lmax + 1)))
#                                      * q2d(np.arange(lmax + 1))
#                                      * wiener)
# dmap_recon = hp.alm2map(dlm_recon, nside=1024)


# # plot
# view_map(dmap_input * mask_b, '', 'dmap_input.png', min=-0.0024, max=0.0024)
# view_map(dmap_recon * mask_b, '', 'dmap_recon.png', min=-0.0024, max=0.0024)























