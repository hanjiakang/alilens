import os
import sys
import healpy as hp

sys.path.insert(0, './')
from one import ALILENS

class simsLensing:
    def __init__(self):
        self.cmbs = os.path.join( ALILENS, 'sims/cmbs/map_TQU_1024_%04d.fits')

    def hashdict(self):
        return {'cmbs': self.cmbs}

    def get_sim_tmap(self, idx):
        noise = hp.read_map( os.path.join(ALILENS, 'sims/noise/map_noise_nside1024_%04d.fits') % idx, field=0)
        return hp.read_map(self.cmbs % idx, field=0) + noise

    def get_sim_pmap(self, idx):
        Q, U = hp.read_map(self.cmbs % idx, field=(1,2))
        noise_Q, noise_U = hp.read_map( os.path.join(ALILENS, 'sims/noise/map_noise_nside1024_%04d.fits') % idx, field=(1,2))
        return Q + noise_Q, U + noise_U
