import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import pymaster as nm

print(1)
mask=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits")
print(2)
aposcale=2
print(3)
apomask = nm.mask_apodization(mask, aposcale, apotype="C1")
print(4)
hp.write_map("apomask.fits",apomask, overwrite=True)

print(1)
mask=hp.read_map("/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits")
print(2)
aposcale=2
print(3)
apomask = nm.mask_apodization(mask, aposcale, apotype="C1")
print(4)
hp.write_map("apomask_b.fits",apomask, overwrite=True)
