from importlib.machinery import SourceFileLoader
import argparse
# import multiprocessing as mp

par = SourceFileLoader('parfile', '/home/jkhan/2021/ali_cpt_lens/params.py').load_module()

parser = argparse.ArgumentParser(description='AliCPT MF calculation')
parser.add_argument('-k', dest='k', action='store', default=[], nargs='+',
                    help='QE keys (NB: both gradient and curl are calculated at the same time)')
args = parser.parse_args()


ivfsA = par.qcls_dd.qeA.f2map1.ivfs
ivfsB = par.qcls_dd.qeB.f2map1.ivfs
ftlA = ivfsA.get_ftl()
felA = ivfsA.get_fel()
fblA = ivfsA.get_fbl()
ftlB = ivfsB.get_ftl()
felB = ivfsB.get_fel()
fblB = ivfsB.get_fbl()


for k in args.k:
    par.n1_dd.get_n1(k, 'p', par.cl_unl['pp'], ftlA, felA, fblA, 2048,
                     kB=k, ftlB=ftlB, felB=felB, fblB=fblB)

# pool = mp.Pool(processes=3)
# for k in args.k:
#     pool.apply_async(par.n1_dd.get_n1, 
#                      args=(k, 'p', par.cl_unl['pp'], ftlA, felA, fblA, 2048,), 
#                      kwds={'kB':k, 'ftlB':ftlB, 'felB':felB, 'fblB':fblB})

# pool.close()
# pool.join()
