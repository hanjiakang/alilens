### 2022 version
This is the 2022 version of AliCPT lensing pipline. For more details, please check out [Ali lensing paper 2022](https://arxiv.org/abs/2204.08158).

To run the code, please first check out *constant.py* and setup your library file path, output directory, dataset information and so on. There is a integrated shell script *one.sh* for this pipeline, which you can run through:

```
  chmod +x one.sh
  
  ./one.sh
```


Then, run *sims.py* to generate a serial of data sets containing lensed CMB signals and inhomogeneous instrumental noises. The lensed CMB signals is calculated using planck FPP potential power spectrum and unlensed TEB CMB signal power spectrum. THe inhomogeneous intrumental noises is generated based on AliCPT noise variance map. A inverse variance filter will be calculated simultaneously.

After that, you can run *plot.py* and start the reconstruction of the lensing potential based on the dataset generated in the previous step.

After the reconstruction process finish, there will be several graphic outputs in the ALILENS directory, among which the  *recon_cl.pdf*  is the power spectrum output and *recon_snr.pdf* is the SNR for each ell bins. Meanwhile, you can run *product.py* and then *klmplot.py* to get a weiner filtered reconstruction potential map as well as the potnetial map input for conparison, I strongly recommand you check this reconstruction potentail map to see if the reconstruction process is carried out correctly.
