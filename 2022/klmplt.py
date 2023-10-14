#Reconstruction map
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
mask=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits")
maske=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/apomask.fits")

def calc_fsky(masks,mask): # fsky calculation
    pixarea=hp.nside2pixarea(1024)
    print(pixarea)
    ret2 = np.ones_like(masks)
    ret4 = np.ones_like(masks)
    ret2 *= masks**2
    ret4 *= masks**4
    order2=np.sum(ret2)
    order4=np.sum(ret4)
    fsky=pixarea/4/np.pi
    fsky=fsky*order2**2/order4
    return fsky

def view_map(m, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
     """ View map.
     """
     # TODO beautify this plot
     rot = [180, 60, 0]


     m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
     m[ m==0. ] = np.nan # in case the input map is an apodization mask

     if min==None: min = m[ ~np.isnan(m) ].min()
     if max==None: max = m[ ~np.isnan(m) ].max()

     hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
     hp.graticule()
     plt.savefig(savePath, dpi=300)


mask_b = np.where(hp.read_map(os.path.join(ALILENS, 'sims/ninv/ninv_t.fits')) > 0, 1, 0)
lmax = 2048
q2k = lambda l: l*(l + 1) / 2 # potential -> convergence
q2d = lambda l: (l*(l + 1)) ** 0.5 # potential -> deflection
cut = np.where((np.arange(lmax + 1) > 8) * (np.arange(lmax + 1) < 2000), 1, 0) # band limit

# wiener filter
wiener_dat = np.loadtxt(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/nlkk.dat')).transpose()
wiener = (wiener_dat[2] - wiener_dat[1]) * utils.cli(wiener_dat[2])

# input deflection map
qlm_input = hp.map2alm(hp.read_map(os.path.join(ALILENS, 'sims/cmbs/map_TQU_%d_%04d.fits'%(2048,nsims))))
dlm_input = hp.almxfl(qlm_input, cut * q2d(np.arange(lmax + 1)))
dmap_input = hp.alm2map(dlm_input, nside=1024)

# reconstruction map
klm_recon = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/dat_klm.fits'))
dlm_recon = hp.almxfl(klm_recon, cut * utils.cli(q2k(np.arange(lmax + 1)))
                                     * q2d(np.arange(lmax + 1))
                                     * wiener)
dmap_recon = hp.alm2map(dlm_recon, nside=1024)

# reconstruction map
klm_mf = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/mf_klm.fits'))
dlm_mf = hp.almxfl(klm_mf, cut * utils.cli(q2k(np.arange(lmax + 1)))
                                     * q2d(np.arange(lmax + 1))
                                     * wiener)
dmap_mf = hp.alm2map(dlm_mf, nside=1024)

# plot
# plot
#view_map(hp.read_map(os.path.join(ALILENS, 'sims/cmbs/map_P_1024_0049.fits'),field=0),'','map_input')
#view_map(hp.read_map(os.path.join(ALILENS, "/disk1/home/hanjk/2022/fg4lens/sims/cmbs/map_TQU_1024_0049.fits"),field=0),'','map_input')
view_map(dmap_input * mask_b, '', 'dmap_input.png', min=-0.0024, max=0.0024)
view_map(dmap_recon * mask_b, '', 'dmap_recon.png', min=-0.0024, max=0.0024)
#view_map(dmap_mf * mask_b, '', 'dmap_mf.png', min=-0.0024, max=0.0024)
#view_map((dmap_recon-dmap_mf) * mask_b, '', 'dmap_dat.png', min=-0.0024, max=0.0024)
