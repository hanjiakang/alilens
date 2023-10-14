#the qe keys you want for the output of lensing reconstruction
lib_qe_keys={
#            'ptt':['royalblue', 'TT'] # Temperature only
           'p_p':['green', 'Pol'] # Polarization only
#           ,'p':['tomato', 'MV']# MV estimator
          } 


ALILENS = "/disk1/home/hanjk/test/2022" #the library directory of your dataset and the final output

nside = 1024
lmax = 2048
dlmax = 1024
seeds = [_ for _ in range(199)]
nlev = [#"/home/jkhan/2021/data/Noise_ALI_IHEP_20200730_48/I_NOISE_95_C_1024.fits",
"/disk1/home/hanjk/Noise_ALI_IHEP_20200730/I_NOISE_150_C_1024.fits"
        ]
fwhm_f = [#19.
         11.
        ]


#library that you store your mask
mask_bb="/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits" #bb mask
mask_apodiz="/disk1/home/hanjk/2022/fg4lens/masks/apomask.fits" # apodized bb mask


#library that you store your 48 module noise data
lib_path_48="/disk1/home/hanjk/2022/fg4lens/sims/noise_48/" #total noise of 48 module
lib_cov_48="/disk1/home/hanjk/2022/fg4lens/sims/conv_48/" #total noise saved for covariance matrix calculation

bias=44 # how many sets of CMB maps used to calculate mean field
var=154 # how many sets of CMB maps used to calculate covariance matrix
nset=9 # the number of gourp
nwidth=22 # the width of each gourp (how many sets in each group) p.s. as the sets in bias wille be to calculate 2 mean field, I recommand that the group width equals to one half of bias.
#library that you store your teb maps as well as your finel recontruction results

#library that you store your teb maps as well as your final recontruction results
lib_cls_res='/disk1/home/hanjk/2022/runs-48/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim%d.fits' #the library directory storing detection noise residue maps
lib_cls_noise='/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%d.fits' 
#the library directory storing total residue maps which is the summation of the detection noise residue maps and foreground residue maps
#lib_cls_con='/disk1/home/hanjk/2022/runs/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%02d.fits'
lib_cls_con='/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%d.fits'
#the noise maps which are used to caculate inverse variance map

#the one set of teb noise map to be used in delensing estimation
lib_Ali_map_noise="/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim198.fits"
lib_Ali_map_res="/disk1/home/hanjk/2022/runs-48/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim198.fits"
