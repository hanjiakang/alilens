# AliCPT Lensing Scripts
This is the 2023 version of AliCPT lensing pipline. To run the code, please first check out *constant.py* and setup your library file path, output directory, dataset information and so on.

Then you can run *transfer.py* to which transfer the pre-existing TEB dataset (CMB maps after foreground removing through ILC method) to TQU maps.

After that, I recommend you to run *sims.py* and *cal_cinv.py*, which will genrate realizaions of CMB signal and calculate inverse variance map respectively.

Both python files above can be accelarated using mp if you add -np (any number) at the end of the command line.

And finally, you can produce lensing reconstruction results using *plot.py* and the powerspectrum will be saved in the directory you put in the *constant.py. There will be *recon_cl.pdf*  which is the power spectrum output and *recon_snr.pdf* which is the SNR for each ell bins.

If you want to test the delensning efficiency of your reconstruction results, there is delensing.py which can be used to estimate the results through CIB internal delensing and CIB-lensing cross delensing.

## SCRIPTS

- *sims.py*: alicpt lensing simulation script use Julien's lenspyx.
        
        python -W ignore sims.py -np 70

  where -np is an argument means number of processors.

- *ali2020_sims.py*: alicpt version of plancklens.planck2018_sims.py used to read observed CMB TQU fits data.
- *transfer.py*: transfer from teb dataset to tqu dataset.
        python transfer.py -np 70
- *cal_cinv.py*: calculate the inverse variance map to be used for inverse-variance filter
- *delensing.py*: delensing estimation pipeline. Out put will be all.pdf for CIB-lensing cross delensing and internal.pdf for internal delensing. The data for these curves in the pdf will be saved in folder ../primordial/
  
