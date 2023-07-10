""" ALL Plot Scripts

    * recon_cl plot
    * SNR plot
    * reconstruction map

"""


import os
import sys
import numpy as np
import healpy as hp
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle # Error boxes
from matplotlib.collections import PatchCollection # Error boxes
from plancklens import utils
from plancklens.qresp import get_response
import lenspyx
sys.path.insert(0, './')
from one import *
import params as par
import bandpowers
import array as arr

# Error boxes
def make_error_boxes(ax, xdata, ydata, xerr, yerr, facecolor, edgecolor='None', alpha=0.5):
    errorboxes = []
    for x, y, xe, ye in zip(xdata, ydata, xerr.T, yerr.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)





# Parameters
savePath = './'
qe_keys = lib_qe_keys
btype = 'agr2'




# Plot Reconstructed power spectrum
fig, ax = plt.subplots()

for key in qe_keys.keys():
    print(key)
    bp = bandpowers.ffp10_binner(key, key, par, btype, ksource='p')
    bells = bp.bin_lavs
#    bpower = (bp.get_dat_bandpowers() - bp.get_rdn0()) * bp.get_bmmc()
#    bpower =np.abs((bp.get_dat_bandpowers() - bp.get_rdn0() - bp.get_n1()) * bp.get_bmmc())
    bpower = (bp.get_dat_bandpowers() - bp.get_mcn0() - bp.get_n1()) * bp.get_bmmc()
    bpower4 = (bp.get_dat_bandpowers() - bp.get_mcn0() - bp.get_n1()) * bp.get_bmmc()
#    bpower = bp.get_n1()
    bpower1 = bp.get_mcn0()
    bpower2 = bp.get_rdn0()
    bpower3 =bp.get_dat_bandpowers()
    bpower5 =bp.get_bmmc()
    
    # Potential Estimator
    yerr = np.sqrt(bp.get_cov().diagonal())
    l, u, c = bandpowers.get_blbubc(btype)
    xerr = np.stack([bp.bin_lavs - l, u - bp.bin_lavs + 1])
    make_error_boxes(ax, bells, bpower, xerr, np.stack([yerr, yerr]), facecolor=qe_keys[key][0], edgecolor='None', alpha=0.2)
    ax.scatter(bells, bpower, label='$C_L^{\hat\phi\hat\phi}-\hatN_L^{(0)}$  (%s)' % qe_keys[key][1], c=qe_keys[key][0], s=13)
    plt.plot(bells,np.zeros(len(bells)),color="k",linestyle="dotted")
 #   ax.scatter(bells, bpower, label='n1', s=13)
 #   ax.scatter(bells, bpower1, label='mcn0', s=13)
 #   ax.scatter(bells, bpower2, label='rdn0', s=13)
 #   ax.scatter(bells, bpower3, label='data', s=13)
 #   ax.scatter(bells, bpower4, label='recondata', s=13)
#    ax.plot(bells, bpower5, label='$bmmc$  (%s)' % qe_keys[key][1],c=qe_keys[key][0])
#    ax.plot(bells, bp.get_rdn0()* bp.get_bmmc(),label='$rdn0*bmmc$  (%s)' % qe_keys[key][1], c=qe_keys[key][0],linestyle='dotted')
    print(np.stack([yerr, yerr]))

# fiducial lensing power spectrum
ax.plot(np.arange(2049), bp.clkk_fid, label=r'$C_L^{\phi\phi, fid}$', c='k')
#ax.semilogx()
#ax.semilogy()
ax.set_title('Lensing Reconstruction', fontsize=12)
ax.set_xlabel(r'$L$', fontsize=9)
ax.set_ylabel(r'$10^7 L^2 (L + 1)^2 C_L^{\phi\phi}/2\pi$', fontsize=9)
ax.legend(loc='upper right', fontsize=7)
#ax.set_xlim(l.min(), u.max())
ax.set_xlim(20, 340)
ax.set_ylim(-1, 3.)
plt.grid()
fig.savefig(os.path.join(ALILENS, 'recon_cl.pdf'), dpi=800)
#fig.savefig('recon_cl.png', dpi=300)

# Plot Reconstructed power spectrum
fig, ax = plt.subplots()

for key in qe_keys.keys():
    print(key)
    bp = bandpowers.ffp10_binner(key, key, par, btype, ksource='p')
    bells = bp.bin_lavs
    bpower = bp.get_dat_bandpowers()  * bp.get_bmmc()
    bpower1= bp.get_rdn0() * bp.get_bmmc()
#    bpower = (bp.get_dat_bandpowers() - bp.get_rdn0())
    
    # Potential Estimator
    yerr = np.sqrt(bp.get_cov().diagonal())
    l, u, c = bandpowers.get_blbubc(btype)
    xerr = np.stack([bp.bin_lavs - l, u - bp.bin_lavs + 1])
    ax.scatter(bells, np.abs(bpower), label='$C_L^{\hat\phi\hat\phi}-\hat N_L^{(0)}$  (%s)' % qe_keys[key][1], c=qe_keys[key][0], s=13)
#    ax.scatter(bells, np.abs(bpower1), label='$C_L^{\hat\phi\hat\phi}-\hat N_L^{(0)}$  (%s)' % qe_keys[key][1], c=qe_keys[key][0], s=13)
    

# fiducial lensing power spectrum
#ax.plot(np.arange(2049), bp.clkk_fid, label=r'$C_L^{\phi\phi, fid}$', c='k')
#ax.semilogx()
#ax.semilogy()
ax.set_title('Lensing Reconstruction', fontsize=12)
ax.set_xlabel(r'$L$', fontsize=9)
ax.set_ylabel(r'$10^7 L^2 (L + 1)^2 C_L^{\phi\phi}/2\pi$', fontsize=9)
ax.legend(loc='lower left', fontsize=7)
ax.set_xlim(l.min(),400)
#ax.set_ylim(0.1, 3.)
fig.savefig(os.path.join(ALILENS, 'reccl.pdf'), dpi=500)
#fig.savefig('recon_cl.png', dpi=300)




# Plot SNR plot
fig, ax = plt.subplots()
for qe_key in qe_keys.keys():
    print(qe_key)
    bp = bandpowers.ffp10_binner(qe_key, qe_key, par, btype, ksource='p')
    cov = bp.get_cov()
    signal = (bp.get_dat_bandpowers() - bp.get_rdn0() - bp.get_n1()) * bp.get_bmmc()
#    signal = (bp.get_dat_bandpowers() - bp.get_rdn0()) * bp.get_bmmc()
#    signal = bp.get_fid_bandpowers()
    l, u, c = bandpowers.get_blbubc(btype)
    SNRs = signal / np.sqrt(cov.diagonal())

    ax.plot(bp.bin_lavs, SNRs, c=qe_keys[qe_key][0], label='SNR (%s)' % qe_keys[qe_key][1])
    print(np.dot(np.dot(signal[1:7], np.matrix(cov[1:7,1:7]).I), signal[1:7])[0,0] ** 0.5)
    print(np.sqrt(sum(signal**2 / cov.diagonal())))
#    print(cov)

ax.semilogx()
ax.set_title('Reconstruction SNR', fontsize=12)
ax.set_xlabel(r'$L$', fontsize=12)
ax.set_ylabel('SNR', fontsize=12)
ax.legend()
ax.set_xlim(20,340)
#

#ax.set_ylim(0.1, 3.)
plt.grid()
fig.savefig(os.path.join(ALILENS, 'recon_snr.pdf'), dpi=500)



# Plot cov plot
fig, ax = plt.subplots()
for qe_key in qe_keys.keys():
    print(qe_key)
    bp = bandpowers.ffp10_binner(qe_key, qe_key, par, btype, ksource='p')
    cov = bp.get_cov()
    bells = bp.bin_lavs
    signal = (bp.get_dat_bandpowers() - bp.get_rdn0() - bp.get_n1()) * bp.get_bmmc()
#    signal = (bp.get_dat_bandpowers() - bp.get_rdn0()) * bp.get_bmmc()
#    signal = bp.get_fid_bandpowers()
    l, u, c = bandpowers.get_blbubc(btype)
    SNRs = signal / np.sqrt(cov.diagonal())
    l, u, c = bandpowers.get_blbubc(btype)
#    print(bells[1:8])
#    multi1=np.tile(bells[1:8]**4,(7,1))
#    multi=multi1*multi1.T
#    print(multi1*multi1.T)
#    matrix= np.log10(multi*cov[1:8,1:8])
    matrix= cov[1:8,1:8].T
#    matrix=np.round_(matrix,decimals=2)
    print(matrix)
    mat=np.tril(matrix)
#    mat=mat-np.diag(mat.diagonal())
#    print(np.diag(mat.diagonal()))
#    mat[mat==0]=False
#    mat=np.round_(mat,decimals=2)
    mask=np.zeros(mat.shape)
    mask[mat==0]=1
    mask[mat!=0]=0
#    mat=array['1.5','','','','','','';
#          '0.16','0.72','','','','','';
#          '-0.1','2.9x10^{-2}','0.34','','','','';
#          '-2.9x10^{-2}','5.2x10^{-2}','0.24','','','';
#          '-2.6x10^{-2}','8.5x10^{-2}','2.3x10^{-2}','0.26','0.26','','';
#          '3.7x10^{-2}','4.4x10^{-2}','2.6x10^{-2}','7.2x10^{-2}','0.32','';
#          '5.6x10^{-2}','7x10^{-2}','5.4x10^{-3}','6.8x10^{-2}','4.3x10^{-2}','0.44';]
 #   print(mask)
 #   matrix=np.log10(matrix-np.min(matrix)+1e13)
 #   print(matrix)   
    sns.set_style('whitegrid')
    fig,ax=plt.subplots()
    ax=sns.heatmap(matrix,annot=mat,mask=mask,cmap=sns.diverging_palette(220,20,n=200),cbar=True)
#    ax=sns.heatmap(matrix,fmt=".2f",cmap=sns.diverging_palette(220,20,n=200))
    ax=sns.heatmap(matrix,fmt=".2f",cmap=sns.diverging_palette(220,20,n=200),cbar=False,xticklabels=bells[1:8].astype(int),yticklabels=bells[8:1].astype(int))
    ax.set_ylabel("central ell of the bin")
    ax.set_xlabel("central ell of the bin")
#    ax.set_xticks(bells[1:7])
#    ax.set_xticklabels(bells[1:7])
#    ax.set_yticks(bells[1:7])
#    ax.set_yticklabels(bells[1:7])
#    plt.setp(ax.get_xticklabels(),ha="center")
#    plt.setp(ax.get_yticklabels(),va="center")
#    im=ax.imshow(matrix,cmap=plt.cm.rainbow)
#    plt.colorbar(im)
    plt.show()
#    print(cov)


    fig.savefig(os.path.join(ALILENS, 'cov.png'), dpi=800)


# Reconstruction map
# def view_map(m, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
#     """ View map.
#     """
#     # TODO beautify this plot
#     rot = [180, 60, 0]


#     m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
#     m[ m==0. ] = np.nan # in case the input map is an apodization mask

#     if min==None: min = m[ ~np.isnan(m) ].min()
#     if max==None: max = m[ ~np.isnan(m) ].max()

#     hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap)
#     hp.graticule()
#     plt.savefig(savePath, dpi=300)


# mask_b = np.where(hp.read_map(os.path.join(ALILENS, 'sims/ninv/ninv_t.fits')) > 0, 1, 0)
# lmax = 2048
# q2k = lambda l: l*(l + 1) / 2 # potential -> convergence
# q2d = lambda l: (l*(l + 1)) ** 0.5 # potential -> deflection
# cut = np.where((np.arange(lmax + 1) > 8) * (np.arange(lmax + 1) < 2000), 1, 0) # band limit

# # wiener filter
# wiener_dat = np.loadtxt(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/nlkk.dat')).transpose()
# wiener = (wiener_dat[2] - wiener_dat[1]) * utils.cli(wiener_dat[2])

# # input deflection map
# qlm_input = hp.map2alm(hp.read_map(os.path.join(ALILENS, 'sims/cmbs/map_P_1024_0300.fits')))
# dlm_input = hp.almxfl(qlm_input, cut * q2d(np.arange(lmax + 1)))
# dmap_input = hp.alm2map(dlm_input, nside=1024)

# # reconstruction map
# klm_recon = hp.read_alm(os.path.join(ALILENS, 'products/COM_Lensing_Inhf_2048_R1/MV/dat_klm.fits'))
# dlm_recon = hp.almxfl(klm_recon, cut * utils.cli(q2k(np.arange(lmax + 1)))
#                                      * q2d(np.arange(lmax + 1))
#                                      * wiener)
# dmap_recon = hp.alm2map(dlm_recon, nside=1024)


# # plot
# view_map(dmap_input * mask_b, '', 'dmap_input.png', min=-0.0024, max=0.0024)
# view_map(dmap_recon * mask_b, '', 'dmap_recon.png', min=-0.0024, max=0.0024)























