import argparse
from importlib.machinery import SourceFileLoader

par = SourceFileLoader('parfile', '/home/jkhan/2021/alicpt_lens/params.py').load_module()

parser = argparse.ArgumentParser(description='AliCPT MF calculation')
parser.add_argument('-k', dest='k', action='store', default=[], nargs='+',
                    help='QE keys (NB: both gradient and curl are calculated at the same time)')
args = parser.parse_args()

for k in args.k:
    par.qlms_dd.get_sim_qlm_mf(k, par.mc_sims_mf_dd[0::2], lmax=par.lmax_qlm)
    par.qlms_dd.get_sim_qlm_mf(k, par.mc_sims_mf_dd[1::2], lmax=par.lmax_qlm)
