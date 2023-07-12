"""To enable complete reconstruction, a parameter file should instantiate

        * the inverse-variance filtered simulation library 'ivfs'
        * the 3 quadratic estimator libraries, 'qlms_dd', 'qlms_ds', 'qlms_ss'.
        * the 3 quadratic estimator power spectra libraries 'qcls_dd, 'qcls_ds', 'qcls_ss'.
          (qcls_ss is required for the MCN0 calculation, qcls_ds and qcls_ss for the RDN0 calculation.
           qcls_dd for the MC-correction, covariance matrix. All three for the point source correction.)
        * the quadratic estimator response library 'qresp_dd'
        * the semi-analytical Gaussian lensing bias library 'nhl_dd'
        * the N1 lensing bias library 'n1_dd'.

    The module bandpowers.py shows how these elements are used to build the reconstructed bandpowers.

    On the first call this module will cache a couple of things will be cached in the directories defined below.

    NOTE Conjugate gradient inversion method is used here instead of Homogeneous filtering or that with rescaling.
         Apodization mask is not a must in Conjugate gradient inversion.
         Plus you dont have to specify the white noise level too.

"""

import os
import healpy as hp
import numpy as np
import sys

import plancklens
from plancklens.filt import filt_util, filt_cinv
from plancklens import utils
from plancklens import qest, qecl, qresp
from plancklens import nhl
from plancklens.n1 import n1
from plancklens.sims import utils as maps_utils
from plancklens.qcinv import cd_solve


sys.path.insert(0, './')
from ali2020_sims import simsLensing
from utils import bl_eft
from one import *


# Data Paths
cls_path = os.path.join(os.path.dirname(os.path.abspath(plancklens.__file__)), 'data', 'cls')
ninv_t_Path = os.path.join(ALILENS, 'sims/ninv/ninv_t.fits')
ninv_p_Path = os.path.join(ALILENS, 'sims/ninv/ninv_p.fits')



# Transfer function and Input power spectrum
cl_unl = utils.camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lenspotentialCls.dat'))
cl_len = utils.camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lensedCls.dat'))


# CMB spectra entering the QE weights (the spectra multplying the inverse-variance filtered maps in the QE legs)
cl_weight = utils.camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lensedCls.dat'))
cl_weight['bb'] *= 0.


# Simulation library for Ali noise
# NOTE in planck there is extra power dcl to better match the data properties
sims = simsLensing()
sims = maps_utils.sim_lib_shuffle(sims, { idx : nsims if idx == -1 else idx for idx in range(-1, nsims) })


# Preconditioner
"""
for [id, pre_ops_descr, lmax, nside, iter_max, eps_min, tr, cache] in chain_descr
by default the chain_descr =

    [[3, ["split(dense(" + pcf + "), 64, diag _cl)"],    256,    128,      3,    0.0, cd_solve.tr_cg, cd_solve.cache_mem()],
     [2, ["split(stage(3),  256, diag_cl)"],            512,    256,      3,    0.0, cd_solve.tr_cg, cd_solve.cache_mem()],
     [1, ["split(stage(2),  512, diag_cl)"],            1024,   512,      3,    0.0, cd_solve.tr_cg, cd_solve.cache_mem()],
     [0, ["split(stage(1), 1024, diag_cl)"],            lmax, nside, np.inf, 1.0e-5, cd_solve.tr_cg, cd_solve.cache_mem()]]

"""
# Julien's suggested eps_min = 1e-3 or 1e-4
chain_descr = [[0, ["diag_cl"], lmax_ivf, nside, np.inf, 1e-3, cd_solve.tr_cg, cd_solve.cache_mem()]]
# Conjuage inversion
ninv_t = [ninv_t_Path]
cinv_t = filt_cinv.cinv_t(libdir_cinvt, lmax_ivf, nside, cl_len, transf, ninv_t, marge_monopole=True, marge_dipole=True, marge_maps=[], chain_descr=chain_descr)
ninv_p = [[ninv_p_Path]]
cinv_p = filt_cinv.cinv_p(libdir_cinvp, lmax_ivf, nside, cl_len, transf, ninv_p, chain_descr=chain_descr)
ivfs_raw = filt_cinv.library_cinv_sepTP(libdir_ivfs, sims, cinv_t, cinv_p, cl_len)

ftl_rs = (np.arange(lmax_ivf + 1) >= lmin_ivf)
fel_rs = (np.arange(lmax_ivf + 1) >= lmin_ivf)
fbl_rs = (np.arange(lmax_ivf + 1) >= lmin_ivf)
ivfs   = filt_util.library_ftl(ivfs_raw, lmax_ivf, ftl_rs, fel_rs, fbl_rs)


# QE libraries instances.
# For the MCN0, RDN0, MC-correction etc calculation, we need in general three of them,
# qlms_dd is the QE library which builds a lensing estimate with the same simulation on both legs
# qlms_ds is the QE library which builds a lensing estimate with a simulation on one leg and the data on the second.
# qlms_ss is the QE library which builds a lensing estimate with a simulation on one leg and another on the seco nd.


# Shuffling dictionary.
# ss_dict remaps idx -> idx + 1 by blocks of 60 up to 300.
ss_dict = { k : v for k, v in zip( np.arange(nsims), np.concatenate([np.roll(range(i*60, (i+1)*60), -1) for i in range(0,5)]))}
ds_dict = { k : -1 for k in range(nsims) }
ivfs_d = filt_util.library_shuffle(ivfs, ds_dict) # always return data map
ivfs_s = filt_util.library_shuffle(ivfs, ss_dict)

qlms_dd = qest.library_sepTP(libdir_qlms_dd, ivfs, ivfs,   cl_len['te'], nside, lmax_qlm=lmax_qlm)
qlms_ds = qest.library_sepTP(libdir_qlms_ds, ivfs, ivfs_d, cl_len['te'], nside, lmax_qlm=lmax_qlm)
qlms_ss = qest.library_sepTP(libdir_qlms_ss, ivfs, ivfs_s, cl_len['te'], nside, lmax_qlm=lmax_qlm)




# qecl libraries instances
mc_sims_bias = np.arange(60)
mc_sims_var  = np.arange(60, 300)

# Only qcls_dd needs a mean-field subtraction.
mc_sims_mf_dd = mc_sims_bias
mc_sims_mf_ds = np.array([])
mc_sims_mf_ss = np.array([])

qcls_dd = qecl.library(libdir_qcls_dd, qlms_dd, qlms_dd, mc_sims_mf_dd)
qcls_ds = qecl.library(libdir_qcls_ds, qlms_ds, qlms_ds, mc_sims_mf_ds)
qcls_ss = qecl.library(libdir_qcls_ss, qlms_ss, qlms_ss, mc_sims_mf_ss)

# Semi-analytical Gaussian lensing bias library:
nhl_dd = nhl.nhl_lib_simple(libdir_nhl_dd, ivfs, cl_weight, lmax_qlm)

# N1 lensing bias library:
n1_dd = n1.library_n1(libdir_n1_dd, cl_len['tt'],cl_len['te'],cl_len['ee'])

# QE response calculation library:
qresp_dd = qresp.resp_lib_simple(libdir_qresp, lmax_ivf, cl_weight, cl_len, {'t': ivfs.get_ftl(), 'e':ivfs.get_fel(), 'b':ivfs.get_fbl()}, lmax_qlm)
