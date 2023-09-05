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
from library_parameter import *

apomask=hp.read_map(mask_apodiz,field=0)
parser = argparse.ArgumentParser(description='Simulation for AliCPT Lensing')
parser.add_argument('-np', dest='np', type=int, default=1, help='Number of processes')
args = parser.parse_args()
print(len(apomask))

    
def lensed_cmbs(cls, maske,nside, savePath, fwhm_f=[], nrms_f=None, lmax=4096, dlmax=1024, facres=-1, seed=0):
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
    print(len(maske))
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
    Tlen = lenspyx.alm2lenmap(tlm_unl, [dlm, None], nside, facres=facres, verbose=False)
    Qlen, Ulen = lenspyx.alm2lenmap_spin([elm_unl, None], [dlm, None], nside, 2, facres=facres, verbose=False)
    tlm_len = hp.map2alm(Tlen, lmax=lmax)
    elm_len, blm_len = hp.map2alm_spin([Qlen, Ulen], 2, lmax=lmax)
    
    # Convolution with transfer function
    Tlen = hp.alm2map(hp.almxfl(tlm_len, transf, inplace=True), nside, verbose=False)
    Qlen, Ulen = hp.alm2map_spin([hp.almxfl(elm_len, transf), hp.almxfl(blm_len, transf)], nside, 2, lmax)

    # Save fits File
    hp.write_map(os.path.join(savePath, fname_TQU % (nside, seed)), [Tlen*maske, Qlen*maske, Ulen*maske], overwrite=True)
    hp.write_map(os.path.join(savePath, fname_P % (nside, seed)), Pmap, overwrite=True)


def noise(ncl, maske,fsky,nside, savePath, seed=0, fwhm_f=None):
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
#    nois=np.load(ncl)
    npix = hp.nside2npix(nside)
    nclinput=ncl/fsky

    fname = 'map_noise_nside%d_%04d.fits'
    mt,mq,mu = hp.synfast(nclinput,nside=1024)
    mt=mt*maske
    mq=mq*maske
    mu=mu*maske
    print(2)
#    noise_150=hp.synfast(nois,nside=nside,pol=True)
#    noise_150=noise_150
#    np.random.seed(seed)
#    if isinstance(nlev, list):
#        if len(nlev) == 1: # single channel noise realization
#            m = np.random.normal(size=(3, npix)) * hp.read_map(nlev[0], verbose=False) \
#                * np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1)
#            print(2)
#        else: # multi-channels
            # NOTE 这里是频段合并最直观的方法，但是失败了
            # assert fwhm_f and fwhm_c
            # nrms_f = [ hp.read_map(path, verbose=False) for path in nlev ]
            # mask = apodize_mask(nrms_f[0])
            # mask_i = np.zeros_like(mask)
            # mask_i[mask!=0] = mask[mask!=0] ** -1
            # m_f = [ np.random.normal(size=(3, npix)) * nrms * \
            #         np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1) for nrms in nrms_f]
            # w_f = wl_f(nrms_f, fwhm_f, fwhm_c, pixwin=True)
            # m = hp.alm2map(sum([ np.array([hp.almxfl(alm, w)
            #             for alm in hp.map2alm(m * mask)])
            #             for m, w in zip(m_f, w_f) ]), nside=nside) * mask_i

            # NOTE 这里是用 effective beam 的方法
#            nrms_f = [ hp.read_map(path, verbose=False) for path in nlev ]
#            nlev_c, transf = bl_eft(nrms_f, fwhm_f, lmax=lmax, pixwin=True, ret_nlev=True)
#            mask = apodize_mask(nrms_f[0])
#            mask_i = np.zeros_like(mask)
#            mask_i[mask!=0] = mask[mask!=0] ** -1
#            m = np.random.normal(size=(3, npix)) * uKamin2uKpix(nlev_c, npix) \
#                * np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1) * mask_i ** 0.5
            
#    else: # white noise
#        m = np.random.normal(size=(3, npix)) * uKamin2uKpix(nlev, npix) \
#            * np.array([1, 2 ** 0.5, 2 ** 0.5]).reshape(3,1)
    print(2)
    hp.write_map(os.path.join(savePath, fname % (nside, seed)), [mt,mq,mu], overwrite=True)


def ninv(masked, savePath, fwhm_f=None):
    """ Noise inveres pixel varance for inhomogeneous filtering.

        * nlev      : it dependes ...

        All AliCPT's detectors are polarized and thus for simplicity
        nlev_Q = nlev_T * sqrt(2)
    """
    assert isinstance(nlev, list)
    nrms= hp.read_map(masked[0], verbose=False)
    fname_ninv_t = 'ninv_t.fits'
    fname_ninv_p = 'ninv_p.fits'
    ninv_t = nrms
    ninv_p = nrms

    # TODO Gaussian smoothing the variance map
#    if len(nlev) == 1:
#        ninv_t = 
#        ninv_p = 
#    else:
        # assert fwhm_f and fwhm_c
        # nside = hp.npix2nside(len(nrms_f[0]))
        # mask_b = np.where(nrms_f[0] != 0, 1, 0)
        # w_f = np.array(wl_f(nrms_f, fwhm_f, fwhm_c, pixwin=True))
        # ninv_t = vmaps2vmap_I([ nrms ** 2 for nrms in nrms_f ], w_f, nside) * mask_b
        # ninv_t[ninv_t!=0] = ninv_t[ninv_t!=0] ** -1
        # ninv_p = ninv_t / 2.

#        assert fwhm_f
#        nlev_c, transf = bl_eft(nrms_f, fwhm_f, lmax=lmax, pixwin=True, ret_nlev=True)
#        npix = len(nrms_f[0])
#        mask = apodize_mask(nrms_f[0])
#        ninv_t = uKamin2uKpix(nlev_c, npix) ** -2 * mask
#        ninv_p = ninv_t / 2.
        

    hp.write_map(os.path.join( savePath, fname_ninv_t ), ninv_t, overwrite=True)
    hp.write_map(os.path.join( savePath, fname_ninv_p ), ninv_p, overwrite=True)

def calc_fsky(masks,mask): # fsky calculation
    ret = np.ones_like(masks)
    ret *= mask**2    
    return sum(ret) / len(masks)

if __name__ == '__main__':
    # MultiProcessing
#    ncll=np.load(ncl)
#    fsky=calc_fsky(apomask,apomask)
#    ninv(masked, savePath_ninv, fwhm_f=fwhm_f)
    pool = mp.Pool(processes=args.np)
    for seed in seeds:
        pool.apply_async(lensed_cmbs, args=(cls_in, apomask,nside, savePath_cmbs),
                kwds={'lmax':lmax, 'dlmax':dlmax, 'fwhm_f':fwhm_f, 'nrms_f':nlev, 'seed':seed})
#        pool.apply_async(noise, args=(ncll,apomask,fsky, nside, savePath_noise),
#                kwds={'fwhm_f':fwhm_f, 'seed':seed})

    pool.close()
    pool.join()

#noise(ncl,nside,savePath_noise)
#    pool.apply_async(noise, args=(nlev, nside, savePath_noise),
#                kwds={'fwhm_f':fwhm_f, 'seed':seeds})
