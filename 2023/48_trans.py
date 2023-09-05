#converting 4 module dataset to 48 module dataset
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import pymaster as nm
import os
import sys
import numpy as np
import healpy as hp
import lenspyx
import multiprocessing as mp
import argparse

sys.path.insert(0, './')
from one import *
from library_parameter import *
from utils import uKamin2uKpix, bl_eft, bl, apodize_mask


parser = argparse.ArgumentParser(description='Simulation for AliCPT Lensing')
parser.add_argument('-np', dest='np', type=int, default=1, help='Number of processes')
args = parser.parse_args()

mask1=hp.read_map(mask_bb)
#mask=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/apomask.fits")
path_48=lib_path_48
cov_48=lib_cov_48

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

#hp.mollview(nb,cmap=plt.cm.coolwarm,min=-1,max=1,coord=["G", "C"])
def transf(cls, nside, savePath, seed=0, fwhm_f=None):
    """ Lensed CMB TQU maps random realiazation.

         * cls      : dict which contiains input power spectra
         * nside    : nside
         * savePath : directory where you store the simulations
         * lmax     : max ell
         * dlmax    : dmax ell
         * fwhm     : Full Width Half Maximum (Beam size)
         * seed     : array which contains random seed
    """
    assert 'OMP_NUM_THREADS' in os.environ.keys(), 'Check your env variable OMP_NUM_THREADS'

    fname_TQU = 'map_TQU_%d_%04d.fits'
    fname_P = 'map_P_%d_%04d.fits'
    fname = 'map_noise_nside%d_%04d.fits'
#    cl=cls%(seed)
    
    cls_noise='/disk1/home/hanjk/2022/fg4lens/sims/residue/map_noise_nside%d_%04d.fits'
    cls_res="/disk1/home/hanjk/2022/fg4lens/sims/noise/map_noise_nside%d_%04d.fits"
    # Lensed maps realization. Correlations between PTE not considered.
    # transfer function NOTE pixwin func is included as well
    nt,ne,nb=hp.read_map(cls_noise%(nside,seed),field=(0,1,2))
    rt,re,rb=hp.read_map(cls_res%(nside,seed),field=(0,1,2))
    
    rtl=rt-nt+nt/np.sqrt(12)
    rel=re-ne+ne/np.sqrt(12)
    rbl=rb-nb+nb/np.sqrt(12)

    ntl=nt/np.sqrt(12)
    nel=nt/np.sqrt(12)
    nbl=nt/np.sqrt(12)
    
    path_48="/disk1/home/hanjk/2022/fg4lens/sims/noise_48/"
    print(seed)
    hp.write_map(os.path.join(savePath, fname % (nside, seed)), [rtl,rel,rbl], overwrite=True)
    

def cal_c(cls, nside, savePath, seed=0, fwhm_f=None):
    """ Lensed CMB TQU maps random realiazation.

         * cls      : dict which contiains input power spectra
         * nside    : nside
         * savePath : directory where you store the simulations
         * lmax     : max ell
         * dlmax    : dmax ell
         * fwhm     : Full Width Half Maximum (Beam size)
         * seed     : array which contains random seed
    """
    assert 'OMP_NUM_THREADS' in os.environ.keys(), 'Check your env variable OMP_NUM_THREADS'

    fname_TQU = 'map_TQU_%d_%04d.fits'
    fname_P = 'map_P_%d_%04d.fits'
    fname = 'map_noise_nside%d_%04d.fits'
#    cl=cls%(seed)
    cls_noise='/disk1/home/hanjk/2022/runs/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim%02d.fits'
    cls_res='/disk1/home/hanjk/2022/runs/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%02d.fits'
    # Lensed maps realization. Correlations between PTE not considered.
    # transfer function NOTE pixwin func is included as well
    t,e,b=hp.read_map(cls_noise % (seed),field=(0,1,2))
    rt,re,rb=hp.read_map(cls_res % (seed),field=(0,1,2))
    
    nt=rt-t+t/np.sqrt(12)
    ne=re-e+e/np.sqrt(12)
    nb=rb-b+b/np.sqrt(12)
    
    nt=nt*mask1
    ne=ne*mask1
    nb=nb*mask1

#    ntl,nel,nbl,nte,neb,ntb=hp.anafast([nt,ne,nb],pol=False)
    ntl1=hp.map2alm(nt,pol=False)
    ntl2=hp.map2alm(ne,pol=False)
    ntl3=hp.map2alm(nb,pol=False)
    #apodization needed
    nt,nq,nu=hp.alm2map([ntl1,ntl2,ntl3],nside=1024,pol=True)
    #checking tqu teb interchangability
    nt=nt*mask1
    nq=nq*mask1
    nu=nu*mask1
    print(seed)
    hp.write_map(os.path.join(savePath, fname % (nside, seed)), [nt,nq,nu], overwrite=True)


#cls_res='/disk1/home/hanjk/lens_sims/maps/TEnilc-Bcilc_proj-noise_11arcmin_sim%d.fits'
#cls_noise='/disk1/home/hanjk/lens_sims/maps/tot-residual_TEnilc-Bcilc_11arcmin_sim%d.fits'
#transf(path_48, nside, path_48, seed=10, fwhm_f=None)
#cal_conv(cov_48, nside,cov_48, seed=5, fwhm_f=None)

if __name__ == '__main__':
    # MultiProcessing
    pool = mp.Pool(processes=args.np)
    for seed in seeds:
        pool.apply_async(transf, args=(path_48, nside, path_48),
                kwds={'fwhm_f':fwhm_f, 'seed':seed})
#        pool.apply_async(cal_c, args=(cov_48, nside, cov_48),
#                kwds={'fwhm_f':fwhm_f, 'seed':seed})


    pool.close()
    pool.join()
