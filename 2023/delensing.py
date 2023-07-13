#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:55:16 2022

@author: wuyi
"""

import os
import matplotlib.pyplot as plt
import lenspyx
from lenspyx.utils import camb_clfile
import healpy as hp, numpy as np
import copy
from constant import *

mask=hp.read_map(mask_bb)
maske=hp.read_map(mask_apodiz)
#mask1=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_20uKcut150_C_1024.fits")
#mask2=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPf_invNvar.fits")

import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))
# make sure the version and path is what you expect
lmax=2000
#Set up a new set of parameters for CAMB
pars = camb.CAMBparams(max_l_tensor=2000)
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

#calculate results for these parameters
results = camb.get_results(pars)

#get dictionary of CAMB power spectra
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
for name in powers: print(name)

#You can calculate spectra for different primordial power spectra without recalculating everything
#for example, let's plot the BB spectra as a function of r
pars.set_for_lmax(4000, lens_potential_accuracy=1)
pars.WantTensors = True
results = camb.get_transfer_functions(pars)

rs = [0.01,0.04,0.08,0.1,0.15,0.2,0.3,0.4]

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
    
def get_blbubc(bin_type):
    if bin_type == 'consext8':
        bins_l = np.array([1,14, 32, 64, 128, 175, 220, 265, 310, 355,499,699])
        bins_u = np.array([13,31,63, 127, 174, 219, 264, 309, 354, 400,500,799])
    elif bin_type == 'agr2':
        bins_l = np.array([8, 21, 40, 66, 101, 145, 199, 264, 339, 426, 526, 638, 763, 902])
        bins_u = np.array([20, 39, 65, 100, 144, 198, 263, 338, 425, 525, 637, 762, 901, 2048])
    elif bin_type == 'xdip':
        bins_l = np.array([8, 264, 902])
        bins_u = np.array([263, 901, 2048])
    elif bin_type == 'pdip':
        bins_l = np.array([8, 101, 426])
        bins_u = np.array([100, 425,  2048])
    elif bin_type == 'lowl':
        bins_l = np.array([2,7])
        bins_u = np.array([8,40])
    elif bin_type == '1_10_unb':
        bins_l = np.arange(1, 11)
        bins_u = bins_l
    elif '_' in bin_type:
        edges = np.int_(bin_type.split('_'))
        bins_l = edges[:-1]
        bins_u = edges[1:] - 1
        bins_u[-1] += 1
    else:
        assert 0, bin_type + ' not implemented'
    return bins_l, bins_u, 0.5 * (bins_l + bins_u)
    
def _get_bil(i, L):  # Bin i window function to be applied to cLpp-like arrays as just described
            ret = (fid_bandpowers[i] / vlpp_den[i]) * vlpp_inv[L] * clkk_fid[L] * kswitch[L]
            ret= cl[L]
            ret *= (L >= bin_lmins[i]) * (L <= bin_lmaxs[i])
            return ret
            
def _get_binnedcl(cl):
      bin_lmins, bin_lmaxs, bin_centers = get_blbubc('consext8')
      nbins = len(bin_centers)
      assert len(cl) > bin_lmaxs[-1], (len(cl), bin_lmaxs[-1])
      ret = np.zeros(nbins)
      for i, (lmin, lmax) in enumerate(zip(bin_lmins, bin_lmaxs)):
            ret[i] = np.average(cl[lmin:lmax + 1])
      return ret

def cli(cl):
    """Pseudo-inverse for positive cl-arrays.
    """
    ret = np.zeros_like(cl)
    ret[np.where(cl > 0)] = 1. / cl[np.where(cl > 0)]
    return ret

lmax = 1024 
dlmax = 512  
nside = 1024 
facres = -1
DGL   = 9.1e9
shot  = 1200.
noise_level0 = 15./np.sqrt(12)   # ¦ÌK-arcmin
pixel_size  = (129600./np.pi/12/nside**2)**0.5*60.
noise_level = noise_level0/pixel_size


cls_path = os.path.join(os.path.dirname(os.path.abspath(lenspyx.__file__)), 'data', 'cls')
cl_len   = camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lensedCls.dat'))
cl_unl   = camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lenspotentialCls.dat'))
cib_spec =  np.load('cib_data.npy')
cib_spec = cib_spec[:,:lmax + dlmax+1]

Kappa = hp.fitsfunc.read_map('./input-output.fits',field=(0,1))
fsky = np.sum(Kappa[1]!=0)/len(Kappa[1])

noise=hp.fitsfunc.read_map(lib_Ali_map_res,field=(0,1,2))
bnoi=hp.anafast(noise[2],lmax=lmax)/fsky
Ali_map_noise=hp.fitsfunc.read_map(lib_Ali_map_noise,field=(0,1,2))
#Ali_map_noise=hp.fitsfunc.read_map("/disk1/home/hanjk/lens_sims/maps/TEnilc-Bcilc_proj-noise_11arcmin_sim198.fits",field=(0,1,2))
#Ali_map_noise=hp.fitsfunc.read_map("/disk1/home/hanjk/lens_sims/maps/tot-residual_TEnilc-Bcilc_11arcmin_sim198.fits",field=(0,1,2))

#vmap = hp.read_map('/Users/wuyi/Desktop/I_NOISE_150_C_1024.fits')

Tmap0, Qmap0, Umap0 = hp.synfast([cl_unl['tt'],cl_unl['ee'],cl_unl['bb'],cl_unl['te']], nside=nside, lmax=lmax + dlmax, new=True)
plm, ilm, inlm = hp.synalm([cl_unl['pp'],cib_spec[2],cib_spec[2]+shot+DGL*cli(cib_spec[0])**3.,cib_spec[3],cib_spec[2],cib_spec[3]], lmax=lmax + dlmax, new=True)

dlm0_cn = hp.almxfl(plm, np.sqrt(np.arange(lmax + 1, dtype=float) * np.arange(1, lmax + 2)))

dlm0 =hp.map2alm(Kappa[0], lmax=lmax + dlmax )#/np.sqrt(fsky)

dlm = hp.map2alm(Kappa[1], lmax=lmax + dlmax )#/np.sqrt(fsky)

Talm0, Ealm0, Balm0 = hp.map2alm([Tmap0, Qmap0, Umap0], lmax=lmax + dlmax)
Qmap_len0_cn, Umap_len0_cn  = lenspyx.alm2lenmap_spin([Ealm0, None], [dlm0_cn, None], nside, 2, facres=facres)
#Qmap_len0_cn, Umap_len0_cn  = lenspyx.alm2lenmap_spin([Ealm0, None], [dlm0, None], nside, 2, facres=facres)
Tmap0_cn, Qmap_len0_cn, Umap_len0_cn = hp.smoothing([Tmap0, Qmap_len0_cn, Umap_len0_cn], fwhm=11./60/180.*np.pi )
Talm0_cn, Ealm_len0_cn, Balm_len0_cn = hp.map2alm([Tmap0_cn, Qmap_len0_cn, Umap_len0_cn], lmax=lmax+dlmax)

#mask = (Ali_map_noise[1]!=0)
Tmap_cn = (hp.alm2map(Talm0_cn, nside=nside) + Ali_map_noise[0]) * mask
Emap_cn = (hp.alm2map(Ealm_len0_cn, nside=nside) + Ali_map_noise[1]) * mask
Bmap_cn = (hp.alm2map(Balm_len0_cn, nside=nside) + Ali_map_noise[2]) * mask

Talm0, Ealm0, Balm0 = hp.map2alm([Tmap0, Qmap0, Umap0], lmax=lmax + dlmax)
Qmap_len0, Umap_len0  = lenspyx.alm2lenmap_spin([Ealm0, None], [dlm0, None], nside, 2, facres=facres)
#Qmap_len0=Qmap_len0_cn
#Umap_len0=Umap_len0_cn
Tmap0, Qmap_len0, Umap_len0 = hp.smoothing([Tmap0, Qmap_len0, Umap_len0], fwhm=11./60/180.*np.pi )
Talm0, Ealm_len0, Balm_len0 = hp.map2alm([Tmap0, Qmap_len0, Umap_len0], lmax=lmax+dlmax)

Tmap = (hp.alm2map(Talm0, nside=nside) + Ali_map_noise[0]) * mask
Emap = (hp.alm2map(Ealm_len0, nside=nside) + Ali_map_noise[1]) * mask
Bmap = (hp.alm2map(Balm_len0, nside=nside) + Ali_map_noise[2]) * mask

#mask = (Ali_map_noise[1]!=0)
Tmap_cn = (hp.alm2map(Talm0_cn, nside=nside) + Ali_map_noise[0]) * mask
Emap_cn = (hp.alm2map(Ealm_len0_cn, nside=nside) + Ali_map_noise[1]) * mask
Bmap_cn = (hp.alm2map(Balm_len0_cn, nside=nside) + Ali_map_noise[2]) * mask



fsky = calc_fsky(mask,mask)
print(fsky)
EE_noise    = hp.anafast(Ali_map_noise[1], lmax=lmax+dlmax)/(np.e**(-np.arange(lmax+dlmax + 1, dtype=float) * np.arange(1, lmax+dlmax+2)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
#fil_Wien    = (cl_unl['ee'][:lmax + dlmax+1]*cli(cl_unl['ee'][:lmax + dlmax+1]+EE_noise))
fil_Wien_cn = (cl_unl['ee'][:lmax + dlmax+1]*cli(cl_unl['ee'][:lmax + dlmax+1]+EE_noise))*(cib_spec[3]*cli(cib_spec[2]+shot+DGL*cli(cib_spec[0])**3.))

Ealm_cn = hp.map2alm(Emap_cn, lmax=lmax+dlmax)#/np.sqrt(fsky)
Ealm = hp.map2alm(Emap, lmax=lmax+dlmax)#/np.sqrt(fsky)

#dlm = hp.almxfl(plm, np.sqrt(np.arange(lmax + 1, dtype=float) * np.arange(1, lmax + 2)))#*fil_Wien[:lmax+1])
dlm_cn = hp.almxfl(inlm, np.sqrt(np.arange(lmax + 1, dtype=float) * np.arange(1, lmax + 2))*fil_Wien_cn[:lmax+1])


elm_cn = hp.almxfl(Ealm_cn, 1/(np.e**(-np.arange(lmax+1, dtype=float) * np.arange(1, lmax+2)*(11/180./60.*np.pi)**2/16/np.log(2))))

elm = hp.almxfl(Ealm, 1/(np.e**(-np.arange(lmax+1, dtype=float) * np.arange(1, lmax+2)*(11/180./60.*np.pi)**2/16/np.log(2))))

Qlen, Ulen  = lenspyx.alm2lenmap_spin([elm, None], [dlm, None], nside, 2, facres=facres)
elm_len, blm_len = hp.map2alm_spin([Qlen, Ulen], 2, lmax=lmax)

Qlen_cn, Ulen_cn  = lenspyx.alm2lenmap_spin([elm_cn, None], [dlm_cn, None], nside, 2, facres=facres)
elm_len_cn, blm_len_cn = hp.map2alm_spin([Qlen_cn, Ulen_cn], 2, lmax=lmax)
ell = np.arange(2, lmax + 1, dtype=int)


Bmap_delens = (Bmap - hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi))*mask
Bmap_delens_cn = (Bmap_cn - hp.smoothing(hp.alm2map(blm_len_cn,nside=nside), fwhm=11./60/180.*np.pi))*mask

Nmap = np.random.normal(0, noise_level, len(Tmap0))

BB_delens=hp.anafast(Bmap_delens,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
BB_lens=1.05*hp.anafast(Bmap)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
np.savez("/disk1/home/hanjk/2022/alicpt_lens/primordial/internal",BB_delens=BB_delens,BB_lens=BB_lens)
#plt.plot(ell, hp.alm2cl(blm_len)[ell], label=r'$\hat{C}_\ell^{\rm BB}$')
#plt.plot(ell, hp.alm2cl(blm_len_cn)[ell]/fsky, label=r'$\hat C_\ell^{\rm BB,cib}$')
#plt.plot(ell, hp.anafast(hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(10.4/180./60.*np.pi)**2/16/np.log(2)))**2, label=r'$C_\ell^{\rm BB,lensed}$')
plt.plot(ell, BB_delens-bnoi[2:], label=r'$C_\ell^{\rm BB,Delensed}$')
plt.plot(ell, BB_lens-bnoi[2:], label=r'$C_\ell^{\rm BB,len}$')
#plt.plot(ell, cl_len['bb'][ell], label=r'$C_\ell^{\rm BB,theory}$')
#plt.plot(ell, hp.anafast(Ali_map_noise[2])[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2, label=r'$Noise$')
#plt.plot(ell, hp.anafast(Nmap)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$N$')
ll=np.arange(lmax+1)+1
for r in rs:
    inflation_params = initialpower.InitialPowerLaw()
    inflation_params.set_params(ns=0.96, r=r)
    results.power_spectra_from_transfer(inflation_params) #warning OK here, not changing scalars
    cl = results.get_tensor_cls(lmax, CMB_unit='muK')
#    print(cl[:,2])
    
    plt.loglog(ll,cl[:,2]/ll/(ll+1),label='r='+str(r))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{BB}\ [\mu K^2]$')

#plt.semilogy()
plt.semilogx()
plt.xlim(13,1000)
plt.ylim(5e-7,2e-5)
plt.legend(fontsize=12,frameon=False)
plt.grid()
plt.savefig('internal.pdf',transparent=True,dpi=300)
plt.close()


BB_delens_cn=hp.anafast(Bmap_delens_cn,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
BB_lens_cn=hp.anafast(Bmap_cn)[ell]/(np.e**(-ell*(ell+1)*(11.7/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
np.savez("/disk1/home/hanjk/2022/alicpt_lens/primordial/cib",BB_delens_cn=BB_delens_cn,BB_lens_cn=BB_lens_cn)
#plt.plot(ell, hp.alm2cl(blm_len_cn)[ell], label=r'$\hat{C}_\ell^{\rm BB}$')
#plt.plot(ell, hp.alm2cl(blm_len)[ell]/fsky, label=r'$\hat C_\ell^{\rm BB,phi}$')
#plt.plot(ell, hp.alm2cl(blm_len_cn)[ell]/fsky, label=r'$\hat C_\ell^{\rm BB,cib}$')
plt.plot(ell, BB_delens_cn-bnoi[2:], label=r'$C_\ell^{\rm BB,Delensed}$')
plt.plot(ell, BB_lens_cn-bnoi[2:], label=r'$C_\ell^{\rm BB,len}$')
#plt.plot(ell, cl_len['bb'][ell], label=r'$C_\ell^{\rm BB,theory}$')
#plt.plot(ell, hp.anafast(Ali_map_noise[2])[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$Noise$')
#plt.plot(ell, hp.anafast(Nmap)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$N$')
ll=np.arange(lmax+1)+1
for r in rs:
    inflation_params = initialpower.InitialPowerLaw()
    inflation_params.set_params(ns=0.96, r=r)
    results.power_spectra_from_transfer(inflation_params) #warning OK here, not changing scalars
    cl = results.get_tensor_cls(lmax, CMB_unit='muK')
#    print(cl[:,2])
    
    plt.loglog(ll,cl[:,2]/ll/(ll+1),label='r='+str(r))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{BB}\ [\mu K^2]$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{BB}\ [\mu K^2]$')

#plt.semilogy()
plt.semilogx()
plt.xlim(13,1000)
plt.legend(fontsize=12,frameon=False)
plt.ylim(5e-7,2e-5)

plt.grid()
plt.savefig('all.pdf',transparent=True,dpi=600)

plt.close()



l=get_blbubc('consext8')[2]
C_cro = hp.anafast(Bmap_cn,hp.smoothing(hp.alm2map(blm_len_cn,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_BB  = hp.anafast(Bmap_cn,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_re  = hp.anafast(hp.smoothing(hp.alm2map(blm_len_cn,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky

plt.plot(l, _get_binnedcl(C_cro**2/C_BB/C_re),label=r'$CIB-cross$')
np.savez("/disk1/home/hanjk/2022/alicpt_lens/primordial/cro",C_cro**2/C_BB/C_re)
C_cro = hp.anafast(Bmap,hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_BB  = hp.anafast(Bmap,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_re  = hp.anafast(hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky

np.savez("/disk1/home/hanjk/2022/alicpt_lens/primordial/in",C_cro**2/C_BB/C_re)
plt.plot(l, _get_binnedcl(C_cro**2/C_BB/C_re),label=r'$Internal$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\rho_{\ell}$')
plt.semilogx()
plt.legend(fontsize=12,frameon=False)
plt.xlim(21,1000)
plt.grid()
plt.savefig("alilens.pdf",transparent=True,dpi=600)
plt.close()

plt.plot(ell, C_cro**2/C_BB/C_re)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\rho_{\ell}$')
plt.semilogx()
plt.savefig("alilens1.png",transparent=True,dpi=300)
plt.close()

C_cro = hp.anafast(Bmap,hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_BB  = hp.anafast(Bmap,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
C_re  = hp.anafast(hp.smoothing(hp.alm2map(blm_len,nside=nside), fwhm=11./60/180.*np.pi)*mask,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky

plt.plot(l, _get_binnedcl(C_cro**2/C_BB/C_re),label=r'$Internal$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\rho$')
plt.semilogx()

plt.savefig("internald.png",transparent=True,dpi=600)
plt.close()

plt.plot(ell, C_cro**2/C_BB/C_re)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\rho$')
plt.semilogx()

plt.savefig("internal1.png",transparent=True,dpi=600)
plt.close()
