""" Utils 0_0
"""
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from plancklens.utils import cli



uKamin2uKpix = lambda n, npix : n / np.sqrt((360 * 60) ** 2 / np.pi / npix)
uKpix2uKamin = lambda n, npix : n * np.sqrt((360 * 60) ** 2 / np.pi / npix)



def bl(fwhm, lmax=None, nside=None, pixwin=True):
    """ Transfer function.

        * fwhm      : beam fwhm in arcmin
        * lmax      : lmax
        * pixwin    : whether include pixwin in beam transfer function
        * nside     : nside
    """
    assert lmax or nside
    lmax = min( 3 * nside - 1, lmax ) if nside and lmax else lmax if lmax else 3*nside - 1
    ret = hp.gauss_beam(fwhm * np.pi / 60. / 180., lmax=lmax)
    if pixwin:
        assert nside is not None
        ret *= hp.pixwin(nside, lmax=lmax)
    return ret


def bl_eft(nrms_f, fwhm_f, lmax=None, pixwin=True, ret_nlev=False):
    """ Effective beam.
    """
    nrms_f = [ hp.read_map(nrms) if isinstance(nrms, str) else nrms for nrms in nrms_f ]
    nside = hp.npix2nside(len(nrms_f[0]))
    nlev_f = np.array([ uKpix2uKamin(np.mean(nrms[nrms > 0] ** -2) ** -0.5 , 
                        hp.nside2npix(nside)) for nrms in nrms_f]) # in uK.radians
    nlev = sum(nlev_f ** -2) ** -0.5 # in uk.arcmin
    bl_f = [ bl(fwhm, pixwin=pixwin, lmax=lmax, nside=nside) for fwhm in fwhm_f ]
    bl_eft = (sum([ nlev ** -2 * bl ** 2 for nlev, bl in zip(nlev_f, bl_f) ])
                * nlev ** 2) ** 0.5
    
    if ret_nlev:
        return nlev, bl_eft
    else:
        return bl_eft
        



def nl(nlev, fwhm, lmax=None, nside=None, pixwin=True):
    """ Detector noise spectrum

        * nlev      : noise level in uK.arcmin
    """
    # uK.arcmin -> uK.radians
    return ( (nlev * np.pi /60. /180.) * \
            cli(bl(fwhm, lmax=lmax, nside=nside, pixwin=pixwin)) ) ** 2


def apodize_mask(m, savePath=None):
    """ Apodization.

        * m         : noise rms map
    """
    # TODO more types of apodizaiton
    mask = np.zeros_like(m)
    mask[m!=0] = m[m!=0] ** -2 * len(m[m>0]) / sum( m[m>0] ** -2 )

    if savePath:
        hp.write_map(savePath, mask, overwrite=True)
        
    return mask

    





def view_map(m, title=None, savePath=None, min=None, max=None, cmap='YlGnBu_r'):
    """ View map.
    """
    # TODO beautify this plot
    rot = [180, 60, 0]


    m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
    # FIXME DONT change the input map

    if min==None: min = m[ ~np.isnan(m) ].min()
    if max==None: max = m[ ~np.isnan(m) ].max()

    hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
    hp.graticule()
    plt.savefig(savePath, dpi=300)



def wl_f(nrms_f, fwhm_f, fwhm_c, pixwin=True, nl_c=False):
    """ Combination weights wl for each freq. channel.
        Substitute of SMICA weights in vmaps2vmap.

        * nrms_f    : noise rms maps
        * fwhm_f    : in arcmin
        * fwhm_c    : in arcmin
        * pixwin    : pixwin
        * nl_c      : if True, return combined noise power spectrum

        NOTE
        没能够成功把不同频段的图结合起来，因为 variance map 还是有一些问题。
    """
    assert len(nrms_f) == len(fwhm_f)
    nrms_f = [ hp.read_map(nrms) if isinstance(nrms, str) else nrms for nrms in nrms_f ]
    nside = hp.npix2nside(len(nrms_f[0]))

    # [95GHz, 150GHz]
    # [11.35047730920254, 17.11698785862552]
    nlev_f = [ uKpix2uKamin(np.mean(nrms[nrms > 0] ** -2) ** -0.5 , hp.nside2npix(nside)) 
               for nrms in nrms_f]
    bl_f = [ bl(fwhm, pixwin=pixwin, nside=nside) for fwhm in fwhm_f ]
    bl_c = bl(fwhm_c, pixwin=pixwin, nside=nside)
    nli_f = [ cli(nl(nlev, fwhm, nside=nside, pixwin=pixwin)) # inverse of nl
              for nlev, fwhm in zip(nlev_f, fwhm_f) ]
    wl_f = np.array([ cli(sum(nli_f)) * cli(bl) * bl_c * nli for bl, nli in zip(bl_f, nli_f) ])

    if nl_c:
        return wl_f, cli(sum(nli_f))
    else:
        return wl_f
