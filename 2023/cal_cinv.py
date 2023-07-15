import sys
import numpy as np
import healpy as hp
import lenspyx
import os 
import matplotlib.pyplot as plt

#sys.path.insert(0, './')
from one import *
mask=hp.read_map(mask_bb)

#def calc_fsky(masks): # fsky calculation
#    ret = np.ones_like(masks[0])
#    for mask in masks:
#        ret *= mask
#    return sum(ret) / len(masks[0])

def calc_fsky(masks,mask): # fsky calculation
    pixarea=hp.nside2pixarea(1024)
    ret2 = np.ones_like(masks)
    ret4 = np.ones_like(masks)
    ret2 *= masks**2
    ret4 *= masks**4
    order2=np.sum(ret2)
    order4=np.sum(ret4)
    fsky=pixarea/4/np.pi
    fsky=fsky*order2**2/order4
    return fsky
    
def cal_fsky(masks,mask): # fsky calculation
    ret2 = np.ones_like(masks)
    ret2 *= masks**2
    order2=np.sum(ret2)
    fsky=order2
    return fsky/len(masks)
fsky=calc_fsky(mask,mask)
fsky1=cal_fsky(mask,mask)




seeds = [_ for _ in range(100)]

#fname = 'map_noise_nside%d_%04d.fits'
fname = 'map_noise_nside1024_%04d.fits'
nt,nq,nu=hp.read_map(os.path.join(savePath_conv, fname % (0)),field=(0,1,2))
#savePath_noi="/disk1/home/hanjk/2022/fg4lens/sims/noise/"
#fname = 'map_noise_nside1024_%04d.fits'
#savePath_res="/disk1/home/hanjk/2022/fg4lens/sims/residue/"
#fname = 'map_noise_nside1024_%04d.fits'
savePath_con="/disk1/home/hanjk/2022/fg4lens/sims/conv/"
fname = 'map_noise_nside1024_%04d.fits'
#savePath_noi="/disk1/home/hanjk/2022/runs/total_residuals/"
#fname = 'tot-residual_TEnilc-Bcilc_11arcmin_sim%02d.fits'
#savePath_res="/disk1/home/hanjk/2022/runs/noise_residuals/"
#fname1 = 'TEnilc-Bcilc_proj-noise_11arcmin_sim%02d.fits'


#

datat=np.zeros([len(nt),len(seeds)])
dataq=np.zeros([len(nt),len(seeds)])
datau=np.zeros([len(nt),len(seeds)])


for seed in seeds:
#    datat[:,seed],dataq[:,seed],datau[:,seed]=hp.read_map(os.path.join(savePath_noi, fname % (seed)),field=(0,1,2))-hp.read_map(os.path.join(savePath_res, fname1 % (seed)),field=(0,1,2))
    datat[:,seed],dataq[:,seed],datau[:,seed]=hp.read_map(os.path.join(savePath_conv, fname % (seed)),field=(0,1,2))
    print(os.path.join(savePath_conv, fname % (seed)))
#    datat[:,seed]=hp.remove_dipole(datat[:,seed])
#    dataq[:,seed]=hp.remove_dipole(dataq[:,seed])
#    datau[:,seed]=hp.remove_dipole(datau[:,seed])

meant=np.mean(datat,axis=1)
#meant=np.transpose(meant)
#datat=(datat-np.transpose(np.tile(meant,(len(seeds),1))))
meane=np.mean(dataq,axis=1)
meanb=np.mean(datau,axis=1)

convt=np.var(datat,axis=1)*mask
convq=np.var(dataq,axis=1)*mask
convu=np.var(datau,axis=1)*mask

stdt=np.std(datat,axis=1)*mask
stdq=np.std(dataq,axis=1)*mask
stdu=np.std(datau,axis=1)*mask

maxt=np.sqrt(stdq.max())
mint=np.sqrt(stdq[np.nonzero(stdq)].min())
print("max",maxt)
print("min",mint)
#binrange=np.logspace(np.log10(1e-12),np.log10(1e5),num=17,base=10)
binrange=np.arange(mint,maxt,(maxt-mint)/20)
print(binrange)
x,bar=np.histogram(np.sqrt(stdq[np.nonzero(stdq)]),bins=binrange)
print(bar[:-1],x)
xlb=np.arange(mint,maxt,(maxt-mint)/20)
plt.bar(xlb[:-1],x,width=(xlb[1]-xlb[0]))
plt.ylabel("counts")
plt.xlabel("STD")
#plt.xticks(xlb,bar[:-1])
plt.savefig("./figure/bar.png")



def view_map(m, title, savePath, min=None, max=None, cmap='YlGnBu_r'):
    """ View map.
    """
    # TODO beautify this plot
    rot = [180, 60, 0]


    m = hp.read_map(m, verbose=False) if isinstance(m, str) else m
    m[ m==0. ] = np.nan # in case the input map is an apodization mask

    if min==None: min = m[ ~np.isnan(m) ].min()
    if max==None: max = m[ ~np.isnan(m) ].max()

    hp.orthview(m, title=title, min=min, max=max, rot=rot, half_sky=True, cmap=cmap,norm="log")
    hp.graticule()
    plt.savefig(savePath, dpi=300)
    print(m)
    m[ np.isnan(m) ] = 0. # in case the input map is an apodization mask
    print(m)

cont=convt
conq=convq
conu=convu

view_map(cont,"convt","./figure/vart.png")
view_map(conq,"convq","./figure/varq.png")
view_map(conu,"convu","./figure/varu.png")



convt[convt!=0] = convt[convt!=0] ** -1
convq[convq!=0] = convq[convq!=0] ** -1
convu[convu!=0] = convu[convu!=0] ** -1
view_map(convt,"convt","./figure/convt.png")
view_map(convq,"convq","./figure/convq.png")
view_map(convu,"convu","./figure/convu.png")

fskyt1=cal_fsky(convt,convt)
fskyq1=cal_fsky(convq,convq)
fskyu1=cal_fsky(convu,convu)
print(fsky1/fskyt1,fsky1/fskyq1,fsky1/fskyu1,fsky1)

#convt=convt*np.sqrt(fsky1/fskyt1)
#convq=convq*np.sqrt(fsky1/fskyq1)
#convu=convu*np.sqrt(fsky1/fskyu1)

fskyt1=cal_fsky(convt,convt)
fskyq1=cal_fsky(convq,convq)
fskyu1=cal_fsky(convu,convu)
print(fskyt1,fskyq1,fskyu1)
#convt=hp.smoothing(convt,0.003199768,pol=False,lmax=lmax)*mask
#convq=hp.smoothing(convq,0.003199768,pol=False,lmax=lmax)*mask
#convu=hp.smoothing(convu,0.003199768,pol=False,lmax=lmax)*mask

view_map(convt,"convt","./figure/cont.png")
view_map(convq,"convq","./figure/conq.png")
view_map(convu,"convu","./figure/conu.png")
#view_map(mask,"mask","./figure/mask.png")

hp.write_map(os.path.join(savePath_ninv, 'ninv_t.fits'), convt, overwrite=True)
hp.write_map(os.path.join(savePath_ninv, 'ninv_p.fits'), convq, overwrite=True)
#view_map(meant,"convt","./figure/cont.png")
#view_map(meane,"conve","./figure/cone.png")
#view_map(meanb,"convb","./figure/conb.png")

#cl='/disk1/home/hanjk/2022/runs/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim02.fits'
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/convt.png")


#cl='/disk1/home/hanjk/2022/runs/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim04.fits'
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/convt1.png")

#cl='/disk1/home/hanjk/2022/runs/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim06.fits'
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/convt2.png")

#cl='/disk1/home/hanjk/2022/runs/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim08.fits'
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/convt3.png")


#cl=os.path.join(savePath_res, fname % (1024, 2))
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/nt.png")
#view_map(ne,"nq","./figure/nq.png")
#view_map(nb,"nu","./figure/nu.png")

#cl=os.path.join(savePath_res, fname % (1024, 4))
#nt,ne,nb=hp.read_map(cl,field=(0,1,2))
#view_map(nt,"nt","./figure/nt1.png")
#view_map(ne,"nq","./figure/nq1.png")
#view_map(nb,"nu","./figure/nu1.png")
