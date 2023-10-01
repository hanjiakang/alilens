import os
import sys
import numpy as np
import healpy as hp
import lenspyx
import multiprocessing as mp
import argparse

sys.path.insert(0, './')
from one import *
from utils import uKamin2uKpix, bl_eft, bl, apodize_mask


parser = argparse.ArgumentParser(description='Simulation for AliCPT Lensing')
parser.add_argument('-np', dest='np', type=int, default=1, help='Number of processes')
args = parser.parse_args()

def lensed_cmbs(cls, nside, savePath, fwhm_f=[], nrms_f=None, lmax=4096, dlmax=1024, facres=-1, seed=0):
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
    
    # Lensed maps realization. Correlations between PTE not considered.
    # transfer function NOTE pixwin func is included as well
    if len(fwhm_f) == 0 or 1:
        fwhm_f = fwhm_f[0] if fwhm_f else 0.
        transf = bl(fwhm_f, nside=nside, lmax=lmax, pixwin=True)
    else:
        assert nrms_f
        transf = bl_eft(nrms_f, fwhm_f, nside=nside, lmax=lmax, pixwin=True)

    np.random.seed(seed)

    # Unlensed TQUP alms
    plm = hp.synalm(cls['pp'], lmax=lmax + dlmax, new=True, verbose=False)
    Pmap = hp.alm2map(plm, nside, verbose=False)
    dlm = hp.almxfl(plm, np.sqrt(np.arange(lmax + 1, dtype=float) * np.arange(1, lmax + 2)))
    tlm_unl, elm_unl, blm_unl = hp.synalm(  [cls['tt'], cls['ee'], cls['bb'], cls['te']], lmax=lmax+dlmax, new=True, verbose=False)
    # NOTE we only consider lensing induced E->B modes
    
    # Lensed TQU maps
 geom_info = ('healpix', {'nside':nside}) 
    Tlen = lenspyx.lensing.alm2lenmap(tlm_unl, dlm,  geometry=geom_info, verbose=False)
    #Tlen = lenspyx.alm2lenmap(tlm_unl, [dlm, None], nside, facres=facres, verbose=False)
    #Qlen, Ulen = lenspyx.alm2lenmap_spin([elm_unl, None], [dlm, None], nside, 2, geometry=geom_info, verbose=False)   #######################这里为啥不能用alm2lenmap
    Qlen, Ulen = lenspyx.lensing.alm2lenmap_spin(elm_unl, dlm, 2, geometry=geom_info, verbose=False)
    tlm_len = hp.map2alm(Tlen, lmax=lmax)
    elm_len, blm_len = hp.map2alm_spin([Qlen, Ulen], 2, lmax=lmax)
    
    # Convolution with transfer function
    Tlen = hp.alm2map(hp.almxfl(tlm_len, transf, inplace=True), nside, verbose=False)
    Qlen, Ulen = hp.alm2map_spin([hp.almxfl(elm_len, transf), hp.almxfl(blm_len, transf)], nside, 2, lmax)

    # Save fits File
    hp.write_map(os.path.join(savePath, fname_TQU % (nside, seed)), [Tlen, Qlen, Ulen], overwrite=True)
    hp.write_map(os.path.join(savePath, fname_P % (nside, seed)), Pmap, overwrite=True)
    print(1)


def noise(nlev, nside, savePath, seed=0, fwhm_f=None):
    """ Noise simulations, we include
            - white noise
            - noise realization according to given noise variance map
            - noise realization and combination of multi-channels

        * nlev      : it depends ...
        * nside     : nside
        * savePath  : directory to save data
        * seed      : random seed
        * fwhm_f    : only for multi-channel combination

        All AliCPT's detectors are polarized and thus for simplicity
        nlev_Q = nlev_T * sqrt(2)
    """
    assert isinstance(nlev, list or float or int)
    npix = hp.nside2npix(nside)
    fname = 'map_noise_nside%d_%04d.fits'

    np.random.seed(seed)
    m = np.random.normal(size=(3, npix)) * hp.read_map(nlev[0], verbose=False) \
                * np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1)

    hp.write_map(os.path.join(savePath, fname % (nside, seed)), m, overwrite=True)

def noises(nlev, nside, savePath, seed=0, fwhm_f=None):
    """ Noise simulations, we include
            - white noise
            - noise realization according to given noise variance map
            - noise realization and combination of multi-channels

        * nlev      : it depends ...
        * nside     : nside
        * savePath  : directory to save data
        * seed      : random seed
        * fwhm_f    : only for multi-channel combination

        All AliCPT's detectors are polarized and thus for simplicity
        nlev_Q = nlev_T * sqrt(2)
    """
    npix = hp.nside2npix(nside)
    fname = 'map_noise_nside%d_%04d.fits'

    m = np.random.normal(size=(3, npix)) * uKamin2uKpix(nlev, npix) \
    * np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1)

    hp.write_map(os.path.join(savePath, fname % (nside, seed)), m, overwrite=True)


def ninv(nlev, savePath, fwhm_f=None):
    """ Noise inveres pixel varance for inhomogeneous filtering.

        * nlev      : it dependes ...

        All AliCPT's detectors are polarized and thus for simplicity
        nlev_Q = nlev_T * sqrt(2)
    """
    assert isinstance(nlev, list)
    nrms_f = [hp.read_map(path, verbose=False) for path in nlev]
    fname_ninv_t = 'ninv_t.fits'
    fname_ninv_p = 'ninv_p.fits'


    # TODO Gaussian smoothing the variance map
    if len(nlev) == 1:
        ninv_t = nrms_f[0] ** 2
        ninv_t[ninv_t!=0] = ninv_t[ninv_t!=0] ** -1
        ninv_p = ninv_t / 2.
    else:
        # assert fwhm_f and fwhm_c
        # nside = hp.npix2nside(len(nrms_f[0]))
        # mask_b = np.where(nrms_f[0] != 0, 1, 0)
        # w_f = np.array(wl_f(nrms_f, fwhm_f, fwhm_c, pixwin=True))
        # ninv_t = vmaps2vmap_I([ nrms ** 2 for nrms in nrms_f ], w_f, nside) * mask_b
        # ninv_t[ninv_t!=0] = ninv_t[ninv_t!=0] ** -1
        # ninv_p = ninv_t / 2.

        assert fwhm_f
        nlev_c, transf = bl_eft(nrms_f, fwhm_f, lmax=lmax, pixwin=True, ret_nlev=True)
        npix = len(nrms_f[0])
        mask = apodize_mask(nrms_f[0])
        ninv_t = uKamin2uKpix(nlev_c, npix) ** -2 * mask
        ninv_p = ninv_t / 2.
        

    hp.write_map(os.path.join( savePath, fname_ninv_t ), ninv_t, overwrite=True)
    hp.write_map(os.path.join( savePath, fname_ninv_p ), ninv_p, overwrite=True)
    
#noises(10., 1024, "/disk1/home/hanjk/Noise_ALI_IHEP_20200730_48/sims/", seed=0, fwhm_f=None)

if __name__ == '__main__':
    # MultiProcessing
    pool = mp.Pool(processes=args.np)
    for seed in seeds:
        pool.apply_async(lensed_cmbs, args=(cls_in, nside, savePath_cmbs),
                kwds={'lmax':lmax, 'dlmax':dlmax, 'fwhm_f':fwhm_f, 'nrms_f':nlev, 'seed':seed})
        pool.apply_async(noise, args=(nlev, nside, savePath_noise),
                kwds={'fwhm_f':fwhm_f, 'seed':seed})

    pool.close()
    pool.join()

    ninv(nlev, savePath_ninv, fwhm_f=fwhm_f)
