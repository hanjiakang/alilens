#the qe keys you want for the output of lensing reconstruction
lib_qe_keys={
#            'ptt':['royalblue', 'TT'] # Temperature only
           'p_p':['green', 'Pol'] # Polarization only
#           ,'p':['tomato', 'MV']# MV estimator
          } 


ALILENS = '/home/jkhan/2021/data/150hz_48data' #the library directory of your dataset and the final output


#library that you store your mask
mask_bb=""/disk1/home/hanjk/2022/fg4lens/masks/AliCPT_UNPfg_filled_C_1024.fits"" #bb mask
mask_apodiz=""/disk1/home/hanjk/2022/fg4lens/masks/apomask.fits"" # apodized bb mask


#library that you store your 48 module noise data
lib_path_48="/disk1/home/hanjk/2022/fg4lens/sims/noise_48/" #total noise of 48 module
lib_cov_48="/disk1/home/hanjk/2022/fg4lens/sims/conv_48/" #total noise saved for covariance matrix calculation

bias=44 # how many sets of CMB maps used to calculate mean field
var=154 # how many sets of CMB maps used to calculate covariance matrix

#library that you store your teb maps as well as your finel recontruction results
lib_cls_res='/disk1/home/hanjk/2022/runs-48/noise_residuals/TEnilc-Bcilc_proj-noise_11arcmin_sim%d.fits'
lib_cls_noise='/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%d.fits'
#lib_cls_con='/disk1/home/hanjk/2022/runs/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%02d.fits'
lib_cls_con='/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim%d.fits'

#the one set of teb noise map to be used in delensing estimation
lib_Ali_map_noise="/disk1/home/hanjk/2022/runs-48/total_residuals/tot-residual_TEnilc-Bcilc_11arcmin_sim198.fits"
