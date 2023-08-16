import os
import sys
from lenspyx.utils import camb_clfile
import lenspyx

sys.path.insert(0, './')
from utils import bl_eft
from library_parameter import *



transf = bl_eft(nlev, fwhm_f, lmax=lmax, pixwin=True)

# input power spectrum
cls_path = os.path.join(os.path.dirname(os.path.abspath(lenspyx.__file__)), 'data', 'cls')
cls_in = camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lenspotentialCls.dat'))

# savePath
savePath_cmbs = os.path.join(ALILENS, 'sims', 'cmbs')
savePath_noise = os.path.join(ALILENS, 'sims', 'noise')
savePath_ninv = os.path.join(ALILENS, 'sims', 'ninv')
if not os.path.exists(savePath_cmbs): os.makedirs(savePath_cmbs)
if not os.path.exists(savePath_noise): os.makedirs(savePath_noise)
if not os.path.exists(savePath_ninv): os.makedirs(savePath_ninv)




# reconstruction config
nsims = len(seeds) - 1
lmax_ivf = 2048
lmin_ivf = 100
lmax_qlm = 2048


# savePath
libdir_cinvt    = os.path.join(ALILENS, 'temp', 'cinv_t')
libdir_cinvp    = os.path.join(ALILENS, 'temp', 'cinv_p')
libdir_ivfs     = os.path.join(ALILENS, 'temp', 'ivfs')
libdir_qlms_dd  = os.path.join(ALILENS, 'temp', 'qlms_dd')
libdir_qlms_ds  = os.path.join(ALILENS, 'temp', 'qlms_ds')
libdir_qlms_ss  = os.path.join(ALILENS, 'temp', 'qlms_ss')
libdir_qcls_dd  = os.path.join(ALILENS, 'temp', 'qcls_dd')
libdir_qcls_ds  = os.path.join(ALILENS, 'temp', 'qcls_ds')
libdir_qcls_ss  = os.path.join(ALILENS, 'temp', 'qcls_ss')
libdir_nhl_dd   = os.path.join(ALILENS, 'temp', 'nhl_dd')
libdir_n1_dd    = os.path.join(ALILENS, 'temp', 'n1_dd')
libdir_qresp    = os.path.join(ALILENS, 'temp', 'qresp')
if not os.path.exists(libdir_cinvt  ): os.makedirs(libdir_cinvt  )
if not os.path.exists(libdir_cinvp  ): os.makedirs(libdir_cinvp  )
if not os.path.exists(libdir_ivfs   ): os.makedirs(libdir_ivfs   )
if not os.path.exists(libdir_qlms_dd): os.makedirs(libdir_qlms_dd)
if not os.path.exists(libdir_qlms_ds): os.makedirs(libdir_qlms_ds)
if not os.path.exists(libdir_qlms_ss): os.makedirs(libdir_qlms_ss)
if not os.path.exists(libdir_qcls_dd): os.makedirs(libdir_qcls_dd)
if not os.path.exists(libdir_qcls_ds): os.makedirs(libdir_qcls_ds)
if not os.path.exists(libdir_qcls_ss): os.makedirs(libdir_qcls_ss)
if not os.path.exists(libdir_nhl_dd ): os.makedirs(libdir_nhl_dd )
if not os.path.exists(libdir_n1_dd  ): os.makedirs(libdir_n1_dd  )
if not os.path.exists(libdir_qresp  ): os.makedirs(libdir_qresp  )
