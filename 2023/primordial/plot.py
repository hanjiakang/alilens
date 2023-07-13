import os
import matplotlib.pyplot as plt
import lenspyx
from lenspyx.utils import camb_clfile
import healpy as hp, numpy as np
import copy

#mask=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits")
#maske=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/apomask.fits")
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




def get_blbubc(bin_type):
    if bin_type == 'consext8':
        bins_l = np.array([1,14, 32, 46,64, 80,105,128, 155,175, 220, 265, 310, 355,499,699])
        bins_u = np.array([13,31,45,63,79, 104,127, 154,174, 219, 264, 309, 354, 400,500,799])
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

lmax=1024
#Set up a new set of parameters for CAMB
pars = camb.CAMBparams(max_l_tensor=2000)
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
ell = np.arange(2, lmax + 1, dtype=int)
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


Kappa = hp.fitsfunc.read_map('../delensing/input-output00.fits',field=(0,1))
fsky = np.sum(Kappa[1]!=0)/len(Kappa[1])


cls_path = os.path.join(os.path.dirname(os.path.abspath(lenspyx.__file__)), 'data', 'cls')
cl_len   = camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lensedCls.dat'))
cl_unl   = camb_clfile(os.path.join(cls_path, 'FFP10_wdipole_lenspotentialCls.dat'))
#cib_spec =  np.load('../cib_data.npy')
#cib_spec = cib_spec[:,:lmax + dlmax+1]


noise=hp.fitsfunc.read_map("/disk1/home/hanjk/2022/runs-48/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim198.fits",field=(0,1,2))
bnoi=hp.anafast(noise[2],lmax=lmax)/fsky

datainternal1=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/internal1.npz")
datacib1=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cib1.npz")

datainternal2=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/internal2.npz")
datacib2=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cib2.npz")

datainternal3=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/internal3.npz")
datacib3=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cib3.npz")

BB_delens1=datainternal1["BB_delens"]
BB_lens1=datainternal1["BB_lens"]

BB_delens_cn1=datacib1["BB_delens_cn"]
BB_lens_cn1=datacib1["BB_lens_cn"]

BB_delens2=datainternal2["BB_delens"]
BB_lens2=datainternal2["BB_lens"]

BB_delens_cn2=datacib2["BB_delens_cn"]
BB_lens_cn2=datacib2["BB_lens_cn"]

BB_delens3=datainternal3["BB_delens"]
BB_lens3=datainternal3["BB_lens"]

BB_delens_cn3=datacib3["BB_delens_cn"]
BB_lens_cn3=datacib3["BB_lens_cn"]



l=get_blbubc('consext8')[2]
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens1-bnoi[2:])/2/np.pi/1.08), color='b',linestyle="-",label=r'$C_\ell^{\rm BB,Delensed}1$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens1-bnoi[2:])/2/np.pi/1.08), color='b',linestyle="--",label=r'$C_\ell^{\rm BB,len}1$')

plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens2-bnoi[2:])/2/np.pi/1.08), color='r',linestyle="-",label=r'$C_\ell^{\rm BB,Delensed}2$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens2-bnoi[2:])/2/np.pi/1.08),color='r',linestyle="--",label=r'$C_\ell^{\rm BB,len}2$')

plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens3-bnoi[2:])/2/np.pi/1.08), color='g',linestyle="-",label=r'$C_\ell^{\rm BB,Delensed}3$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens3-bnoi[2:])/2/np.pi/1.08), color='g',linestyle="--",label=r'$C_\ell^{\rm BB,len}3$')
#plt.plot(ell, ell*(ell+1)*cl_len['bb'][ell]/2/np.pi, label=r'$C_\ell^{\rm BB,theory}$')
#plt.plot(ell, ell*(ell+1)*bnoi[2:]/2/np/pi)
#plt.plot(ell, hp.anafast(Ali_map_noise[2])[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2, label=r'$Noise$')
#plt.plot(ell, hp.anafast(Nmap)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$N$')


ll=np.arange(lmax+1)+1
rs = [0.032]
for r in rs:
    inflation_params = initialpower.InitialPowerLaw()
    inflation_params.set_params(ns=0.96, r=r)
    results.power_spectra_from_transfer(inflation_params) #warning OK here, not changing scalars
    cl = results.get_tensor_cls(lmax, CMB_unit='muK')
#    print(cl[:,2])
    
    plt.loglog(ll,cl[:,2],color="k",label='r='+str(r))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{BB}/ (2\pi \mu{\rm K}^2)$')

cl1 = results.get_lensed_scalar_cls(lmax, CMB_unit='muK')
#plt.semilogy()
plt.semilogx()
plt.xlim(13,300)
plt.ylim(3e-4,1e-2)
plt.legend(fontsize=12,frameon=False)
#plt.loglog(ll,cl1[:,2],label='r='+str(r))
plt.grid(which='both')
plt.savefig('internal.pdf',transparent=True,dpi=300)
plt.close()


#BB_delens_cn=hp.anafast(Bmap_delens_cn,lmax=lmax+dlmax)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
#BB_lens_cn=hp.anafast(Bmap_cn)[ell]/(np.e**(-ell*(ell+1)*(11.7/180./60.*np.pi)**2/16/np.log(2)))**2/fsky
l=get_blbubc('consext8')[2]
#plt.plot(ell, hp.alm2cl(blm_len_cn)[ell], label=r'$\hat{C}_\ell^{\rm BB}$')
#plt.plot(ell, hp.alm2cl(blm_len)[ell]/fsky, label=r'$\hat C_\ell^{\rm BB,phi}$')
#plt.plot(ell, hp.alm2cl(blm_len_cn)[ell]/fsky, label=r'$\hat C_\ell^{\rm BB,cib}$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens_cn1-1.1*bnoi[2:])/2/np.pi),color='b',linestyle='-', label=r'$C_\ell^{\rm BB,Delensed}$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens_cn1-bnoi[2:])/2/np.pi),color='b', linestyle='--',label=r'$C_\ell^{\rm BB,len}$')

plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens_cn2-1.1*bnoi[2:])/2/np.pi),color='r',linestyle='-', label=r'$C_\ell^{\rm BB,Delensed}$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens_cn2-bnoi[2:])/2/np.pi),color='r', linestyle='--',label=r'$C_\ell^{\rm BB,len}$')

plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_delens_cn3-1.1*bnoi[2:])/2/np.pi),color='g',linestyle='-', label=r'$C_\ell^{\rm BB,Delensed}$')
plt.plot(l, _get_binnedcl(ell*(ell+1)*(BB_lens_cn3-bnoi[2:])/2/np.pi), color='g', linestyle='--',label=r'$C_\ell^{\rm BB,len}$')

#plt.plot(ell, cl_len['bb'][ell], label=r'$C_\ell^{\rm BB,theory}$')

#plt.plot(ell, hp.anafast(Ali_map_noise[2])[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$Noise$')
#plt.plot(ell, hp.anafast(Nmap)[ell]/(np.e**(-ell*(ell+1)*(11/180./60.*np.pi)**2/16/np.log(2)))**2/fsky, label=r'$N$')
ll=np.arange(lmax+1)+1
rs = [0.032,0.02]
col=['--','-']
for i,r in enumerate(rs):
    inflation_params = initialpower.InitialPowerLaw()
    inflation_params.set_params(ns=0.96, r=r)
    results.power_spectra_from_transfer(inflation_params) #warning OK here, not changing scalars
    cl = results.get_tensor_cls(lmax, CMB_unit='muK')
#    print(cl[:,2])
    
    plt.loglog(ll,cl[:,2],color="k",linestyle=col[int(i)],label='r='+str(r))
plt.xlabel(r'$\ell$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{BB}/ (2\pi \mu{\rm K}^2)$')

#plt.semilogy()
plt.semilogx()
plt.xlim(13,300)
plt.legend(fontsize=12,frameon=False)
plt.ylim(3e-4,1e-2)

plt.grid(which='both')
plt.savefig('all.pdf',transparent=True,dpi=300)

plt.close()


C1=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cro1.npz")['arr_0']
C2=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cro2.npz")['arr_0']
C3=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/cro3.npz")['arr_0']
#print(sorted(C1.files))
#plt.plot(l, C1,label=r'$CIB-cross1$')
#plt.plot(l, C2,label=r'$CIB-cross2$')
#plt.plot(l, C3,label=r'$CIB-cross3$')

plt.plot(l, _get_binnedcl(C1),label=r'$CIB-cross1$')
plt.plot(l, _get_binnedcl(C2),label=r'$CIB-cross2$')
plt.plot(l, _get_binnedcl(C3),label=r'$CIB-cross3$')

I1=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/in1.npz")['arr_0']
I2=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/in2.npz")['arr_0']
I3=np.load("/disk1/home/hanjk/2022/alicpt_lens/primordial/in3.npz")['arr_0']
plt.plot(l, _get_binnedcl(I1),label=r'$Internal1$')
plt.plot(l, _get_binnedcl(I2),label=r'$Internal2$')
plt.plot(l, _get_binnedcl(I3),label=r'$Internal3$')

plt.xlabel(r'$\ell$')
plt.ylabel(r'$\rho^2_{\ell}$')
plt.semilogx()
plt.legend(fontsize=12,frameon=False)
plt.xlim(21,1000)
plt.grid(which='both')
plt.savefig('alilens.pdf',transparent=True,dpi=300)

plt.close()


