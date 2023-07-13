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
mask1=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_20uKcut150_C_1024.fits")
mask2=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPf_invNvar.fits")

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

def view_map(m,m1, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
     """ View map.
     """
     #     # TODO beautify this plot
     rot = [180, 60, 0]


     m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
     m[m==0]=hp.UNSEEN
     m1[m1==0]=hp.UNSEEN

     if min==None: min = m[ ~np.isnan(m) ].min()
     if max==None: max = m[ ~np.isnan(m) ].max()
     xsize=1800
     ysize=1400
#     hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.cartview(m, title=title, min=min, max=max, rot=rot, cmap=cmap)
     hp.azeqview(m,xsize=xsize,ysize=ysize,reso=hp.nside2resol(1024,arcmin=True),sub=(1,1,1),rot=rot,min=-0.0024, max=0.0024,cmap=cmap,badcolor='white',title='',return_projected_map=False,cbar=1)
     hp.graticule(color="silver",ls="-",alpha=0.7)
     plt.savefig("dmap.png",dpi=800,bbox_inches='tight')
     plt.figure()
     hp.azeqview(m1,xsize=xsize,ysize=ysize,reso=hp.nside2resol(1024,arcmin=True),sub=(1,1,1),rot=rot,min=min, max=max,cmap=cmap,badcolor='white',title='',return_projected_map=0,cbar=1)
#     hp.azeqview(m1,rot=rot,lamb=True, fig=plt.gcf().number,sub=(1,2,2),cbar=False,cmap='viridis_r',badcolor='white',aspect='auto',title='',return_projected_map=False)
#     hp.orthview(m1, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.cartview(m, title=title, min=min, max=max, rot=rot, cmap=cmap)
     hp.graticule(color="silver",ls="-",alpha=0.7)
     plt.savefig(savePath,dpi=800,bbox_inches='tight')



mask_b = np.where(hp.read_map(os.path.join(ALILENS, 'sims/ninv/ninv_p1.fits')) > 0, 1, 0)
lmax = 2048
q2k = lambda l: l*(l + 1) / 2 # potential -> convergence
q2d = lambda l: (l*(l + 1)) ** 0.5 # potential -> deflection
cut = np.where((np.arange(lmax + 1) > 8) * (np.arange(lmax + 1) < 1500), 1, 0) # band limit

# wiener filter
wiener_dat = np.loadtxt(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/nlkk.dat')).transpose()
wiener = (wiener_dat[2] - wiener_dat[1]) * utils.cli(wiener_dat[2])

# input deflection map
qlm_input = hp.map2alm(hp.read_map(os.path.join(ALILENS, "sims/cmbs/map_P_1024_0198.fits")))
#dlm_input = hp.almxfl(qlm_input, cut) #* q2d(np.arange(lmax + 1)))
dlm_input = hp.almxfl(qlm_input, cut* q2d(np.arange(lmax + 1)))
dmap_input = hp.alm2map(dlm_input, nside=1024)

# reconstruction map
klm_recon = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/dat_klm.fits'))
dlm_recon = hp.almxfl(klm_recon, cut * utils.cli(q2k(np.arange(lmax + 1)))
                                     * q2d(np.arange(lmax + 1))
                                     * wiener
                                     )
dmap_recon = hp.alm2map(dlm_recon, nside=1024)

# reconstruction map
klm_mf = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/mf_klm.fits'))
dlm_mf = hp.almxfl(klm_mf, cut * utils.cli(q2k(np.arange(lmax + 1)))
                                    * q2d(np.arange(lmax + 1))
                                     * wiener
                                     )
dmap_mf = hp.alm2map(dlm_mf, nside=1024)

hp.write_map("input-output.fits",[dmap_input* mask_b*mask,dmap_recon* mask_b*mask],overwrite=True)
#dmap_input,dmap_recon=hp.read_map("input-output.fits",field=(0,1))
#plt.figure()
#fsky1=calc_fsky(maske,maske)
#fsky=calc_fsky(mask,mask)
#cl=hp.anafast([dmap_recon,dmap_input],pol=0)/fsky
#el=np.arange(len(cl[0,:]))
#plt.plot(el,el*(el+1)*cl[0,:]*1e7,label="output")
#plt.plot(el,el*(el+1)*cl[1,:]*1e7,label="input")
#plt.plot(el,el*(el+1)*cl[2,:]*1e7,label="cross")
#plt.xscale("log")
#plt.yscale("log")
#plt.xlabel("ell")
#plt.ylabel("cl")
#plt.xlim(10,2048)
#plt.ylim(1e-2,1e2)
#plt.legend()
#plt.title("cross ps")
#plt.savefig("cross_ps.png")
#plt.legend()
#plt.figure()
#cl=hp.alm2cl(klm_mf)/fsky
#plt.plot(np.arange(len(cl)),cl)
#plt.xscale("log")
#plt.yscale("log")
#plt.xlabel("ell")
#plt.ylabel("cl")
#plt.title("klm2cl_mf")
#plt.savefig("klm2cl.png")
# plot
# plot
#view_map(hp.read_map(os.path.join(ALILENS, 'sims/cmbs/map_P_1024_0049.fits'),field=0),'','map_')
#view_map(hp.read_map(os.path.join(ALILENS, "/disk1/home/hanjk/2022/fg4lens/sims/cmbs/map_TQU_1024_0049.fits"),field=0),'','map_input')
print(dmap_input[dmap_input!=0],dmap_recon[dmap_recon!=0])
plt.figure()
#view_map(dmap_input,dmap_recon, '', 'dmap_input.png', min=-0.0024, max=0.003)

#cmap='binary'
#def view_map(m,m1, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
#     """ View map.
#     """
#     # TODO beautify this plot
#     rot = [180, 60, 0]


#     m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
#     m[ m==0. ] = np.nan # in case the input map is an apodization mask
#     m[m==0]=hp.UNSEEN
#     m1[m1==0]=hp.UNSEEN

#     if min==None: min = m[ ~np.isnan(m) ].min()
#     if max==None: max = m[ ~np.isnan(m) ].max()
#     xsize=1900
#     ysize=1600
#     hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.cartview(m, title=title, min=min, max=max, rot=rot, cmap=cmap)
#     plt.figure()
#     hp.azeqview(m,xsize=xsize,ysize=ysize,reso=hp.nside2resol(1024,arcmin=True),sub=(1,1,1),rot=rot,min=0, max=2,cmap='Greys',badcolor='white',title='',return_projected_map=False,cbar=0)
#     hp.graticule(color="silver",ls="-",alpha=0.7)
#     plt.savefig("mask.png",dpi=800,bbox_inches='tight')
#     plt.figure()
#     fig=hp.azeqview(m1,xsize=xsize,ysize=ysize,reso=hp.nside2resol(1024,arcmin=True),sub=(1,1,1),rot=rot,min=0, max=2,cmap='Greys',badcolor='white',title='',return_projected_map=0,cbar=0)
#     hp.azeqview(m1,rot=rot,lamb=True, fig=plt.gcf().number,sub=(1,2,2),cbar=False,cmap='viridis_r',badcolor='white',aspect='auto',title='',return_projected_map=False)
#    hp.orthview(m1, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.cartview(m, title=title, min=min, max=max, rot=rot, cmap=cmap)
#     hp.graticule(color="silver",ls="-",alpha=0.7)
#     plt.savefig(savePath,dpi=800,bbox_inches='tight')

#view_map(mask,mask2,'','maks2.png', min=0, max=0.07)
#view_map(mask,mask1, '', 'mask1.png', min=0, max=1)
view_map(dmap_recon * mask_b*mask,dmap_input * mask_b*mask, '', 'dmap_recon.png', min=-0.0024, max=0.0024)
#view_map(dmap_input * mask_b*mask, '', 'dmap_input.png', min=-0.0024, max=0.0024)
#view_map(dmap_mf * mask_b*mask, '', 'dmap_mf.png', min=-0.0024, max=0.0024)
#view_map((dmap_recon-dmap_mf) * mask_b*mask, '', 'dmap_dat.png')#, min=-0.0024, max=0.0024)