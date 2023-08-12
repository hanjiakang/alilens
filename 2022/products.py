""" Lensing Products

    Planck 2018 Lensing Products: 
    https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Lensing
    
        COM_Lensing_Inhf_2048_R1
    +       MV/dat_klm.fits
    +       MV/mf_klm.fits
    +       MV/nlkk.dat

        COM_Lensing-SimMap_Inhf_2048_R1
    +       inputs/FFP10_wdipole_lenspotentialCls.dat
    +       MV/sim_klm_{???}.fits
    +       MV/dat_klm.fits
    -       MV/nlkk.dat

"""
import os
import sys
import shutil
import numpy as np
import healpy as hp
import plancklens
from plancklens import utils
import matplotlib.pyplot as plt
from importlib.machinery import SourceFileLoader

sys.path.insert(0, './')
from one import *

par = SourceFileLoader('parfile', "/disk1/home/hanjk/alicpt_lens/params.py").load_module()

PATH = os.path.join(ALILENS, 'products')
if not os.path.exists(PATH): os.makedirs(PATH)



def klm(k, idx, qlms_dd, qresp_dd, savePath, mc_sims_mf=None):
    """ klm products
        qlm (potential) -> klm (convergence)
            klm = l(l+1)/2 qlm

        * i         : ith qlm, if i==-1 return data qlm else sim qlm
        * qlms_dd   : qest library
        * savePath  : savePath
        * mc_sims_mf: idx to calc mf

    """
    lmax = qlms_dd.get_lmax_qlm(k)
    q2k = lambda l : l * (l + 1) / 2
    qresp = qresp_dd.get_response(k, 'p')
    assert len(qresp) == lmax + 1
    klm = hp.almxfl(qlms_dd.get_sim_qlm(k, idx), utils.cli(qresp) * q2k(np.arange(lmax + 1)))

    if mc_sims_mf is not None:
        assert idx==-1
        klm_mf = hp.almxfl(qlms_dd.get_sim_qlm_mf(k, mc_sims_mf), utils.cli(qresp) * q2k(np.arange(lmax + 1)))
        hp.write_alm(os.path.join(savePath, 'dat_klm.fits'), (klm - klm_mf), overwrite=True)
        hp.write_alm(os.path.join(savePath, 'mf_klm.fits'), klm_mf, overwrite=True)
        
    else:
        hp.write_alm(os.path.join(savePath, 'dat_klm.fits' if idx==-1 else 'sim_klm_%03d.fits' % idx), 
                klm, overwrite=True)



def nlkk(k, clpp, nhl_dd, qresp_dd, savePath):
    """ Create nlkk.dat for wiener filtering
        save filter into a text file with cols (L, NL, CL + NL)
        nl is calculated using semi-analytical gaussian lensing bias module `nhl`.

        * k         : key
        * clpp      : potential power spectrum (no weighting)
        * nhl_dd    : nhl library
        * qresp_dd  : resp library
    """
    # TODO we need nlkk.dat for COM_Lensing-SimMap_Inhf like Planck
    lmax = min([len(clpp) - 1, nhl_dd.lmax_qlm])
    q2k = lambda l: l*(l + 1) / 2

    ll = np.arange(lmax + 1)
    nlkk = nhl_dd.get_sim_nhl(-1, k, k)*q2k(ll)** 2*utils.cli(qresp_dd.get_response(k, 'p') ** 2)

    clkk = clpp[:lmax + 1] * q2k(ll) ** 2
    
    # TODO 数据靠右对齐
    np.savetxt(os.path.join(savePath, 'nlkk.dat'), 
            np.array([ll, nlkk, clkk + nlkk]).transpose(),
            fmt = ['%d', '%.12e', '%.12e'],
            header = 'cols = (L, NL, CL + NL)'
            )

def mask():
    pass
    


# COM_Lensing_Inhf
PATH_COM_Lensing = os.path.join(PATH, 'COM_Lensing_Inhf_%d_R1' % lmax)
if not os.path.exists(PATH_COM_Lensing):os.makedirs(PATH_COM_Lensing)

k = 'ptt'
if not os.path.exists(os.path.join(PATH_COM_Lensing, 'MV')):
    os.makedirs(os.path.join(PATH_COM_Lensing, 'MV'))
klm(k, -1, par.qlms_dd, par.qresp_dd, os.path.join(PATH_COM_Lensing, 'MV'), 
        mc_sims_mf=par.mc_sims_mf_dd)
nlkk(k, par.cl_unl['pp'], par.nhl_dd, par.qresp_dd, 
        os.path.join(PATH_COM_Lensing, 'MV'))




# COM_Lensing-SimMap_Inhf
# PATH_COM_Lensing_SimMap = os.path.join(PATH, 'COM_Lensing-SimMap_Inhf_%d_R1' % lmax)
# os.makedirs(PATH_COM_Lensing_SimMap)

# k = 'p'
# os.makedirs(os.path.join(PATH_COM_Lensing_SimMap, 'inputs'))
# os.makedirs(os.path.join(PATH_COM_Lensing_SimMap, 'MV'))

# cls_unl_path = os.path.join(os.path.dirname(os.path.abspath(plancklens.__file__)), 
#         'data', 'cls', 'FFP10_wdipole_lenspotentialCls.dat')
# shutil.copy(cls_unl_path, os.path.join(PATH_COM_Lensing_SimMap, 'inputs'))

# for idx in range(-1, 300):
#     klm(k, idx, par.qlms_dd, par.qresp_dd,
#             os.path.join(PATH_COM_Lensing_SimMap, 'MV'))

