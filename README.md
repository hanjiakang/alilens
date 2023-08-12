<!---
STYLE CONVENTION USED
    bolt italic:
        ***file***"
    code:
       `program` or `library``
       `commands` or `paths`
       `variable`
    bold code:
        **`function`**
        **`type`** or **`structure`**
-->
# alilens

Python 3.9.7

alilens is a python script for cosmology containning most of AliCPT CMB lensing pipeline, developed by Bin Hu, Jiakang Han, Jinyi Liu on behalf of the *AliCPT* collaboration ([publication here.](https://arxiv.org/abs/1710.03047)) This script is mostly based on [Plancklens](https://github.com/carronj/plancklens) by Julien Carron. 

This script may be used for:
* lensing reconstruction on *AliCPT* dataset consists of foregorund residue, CMB signal and instrumental noise
* delensing based on lensing reconstruction and its cross-correlation with Cosmic Microwave Background

### Requirements
To use alilens, you will need [Plancklens](https://github.com/carronj/plancklens), [lenspyx](https://github.com/carronj/lenspyx) and [healpy](https://github.com/carronj/lenspyx) installed locally. Check out their github for instructions. 

### Contents

This code contains most of the AliCPT lensing reconstruction pipeline (see Ali lensing paper [2022](https://arxiv.org/abs/2204.08158), [2023](https://arxiv.org/abs/2303.05705).) Here we provide two versions corresponding to the above papers respectively:

* 2022 version consisting of data simulation scripts, which generate skymaps with planck FPP CMB signal and instrumental noise based on AliCPT noise variance map. It calculates the lensing reconstruciton results for both 95 GHz and 150 GHz from temperature only, polarization only and MV estimators.
* 2023 version consistiing of lensing reconstruction pipeline for polarization only estimation of pre-generated data base, which contains maps with CMB signal, instrumental noise and foreground residues after ILC fg removing method.

### 2022 version
This is the 2022 version of AliCPT lensing pipline. For more details, please check out [Ali lensing paper 2022](https://arxiv.org/abs/2204.08158).

To run the code, please first check out *constant.py* and setup your library file path, output directory, dataset information and so on. There is a integrated shell script *one.sh* for this pipeline, which you can run through:

```
  chmod +x one.sh
  
  ./one.sh
```


Also, you can use seperated python files. Run *sims.py* to generate a serial of data sets containing lensed CMB signals and inhomogeneous instrumental noises. The lensed CMB signals is calculated using planck FPP potential power spectrum and unlensed TEB CMB signal power spectrum. THe inhomogeneous intrumental noises is generated based on AliCPT noise variance map. A inverse variance filter will be calculated simultaneously.

After that, you can run *plot.py* and start the reconstruction of the lensing potential based on the dataset generated in the previous step.

After the reconstruction process finish, there will be several graphic outputs in the ALILENS directory, among which the  *recon_cl.pdf*  is the power spectrum output and *recon_snr.pdf* is the SNR for each ell bins. Meanwhile, you can run *product.py* and then *klmplot.py* to get a weiner filtered reconstruction potential map as well as the potnetial map input for conparison, I strongly recommand you check this reconstruction potentail map to see if the reconstruction process is carried out correctly.

#### SCRIPTS

- *sims.py*: alicpt lensing simulation script use Julien's lenspyx.
        
        python -W ignore sims.py -np 70

  where -np is an argument means number of processors.

- *ali2020_sims.py*: alicpt version of plancklens.planck2018_sims.py used to read observed CMB TQU fits data.

- *klmplot.py*: delensing estimation pipeline. Out put will be all.pdf for CIB-lensing cross delensing and internal.pdf for internal delensing. The data for these curves in the pdf will be saved in folder ../primordial/.

- *product.py*: generate reconstruction klm and semi-analytical n0 nlkk.
### 2023 version
This is the 2023 version of AliCPT lensing pipline. For more details, please check out [Ali lensing paper 2023](https://arxiv.org/abs/2303.05705).

To run the code, please first check out *constant.py* and setup your library file path, output directory, dataset information and so on. There is a integrated shell script *one.sh* for this pipeline, which you can run through:

```
  chmod +x one.sh
  
  ./one.sh
```

Also, you can use seperated python files. You can run *transfer.py* to which transfer the pre-existing TEB dataset (CMB maps after foreground removing through ILC method) to TQU maps.

After that, I recommend you to run *sims.py* and *cal_cinv.py*, which will genrate realizaions of CMB signal and calculate inverse variance map respectively.

*transfer.py* and *sims.py* can be accelarated using mp if you add -np (any number) at the end of the command line.

And finally, you can produce lensing reconstruction results using *plot.py* and the powerspectrum will be saved in the directory you put in the *constant.py. There will be *recon_cl.pdf*  which is the power spectrum output and *recon_snr.pdf* which is the SNR for each ell bins.

If you want to test the delensning efficiency of your reconstruction results, there is *delensing.py* which can be used to estimate the results through CIB internal delensing and CIB-lensing cross delensing. You must provide a reconstruction potential map befroe you run the delensing estimation, which can be caculated using *klmplot.py*.

#### SCRIPTS

- *sims.py*: alicpt lensing simulation script use Julien's lenspyx.
        
        python -W ignore sims.py -np 70

  where -np is an argument means number of processors.

- *ali2020_sims.py*: alicpt version of plancklens.planck2018_sims.py used to read observed CMB TQU fits data.
  
- *transfer.py*: transfer from teb dataset to tqu dataset.
  
        python transfer.py -np 70

- *cal_cinv.py*: calculate the inverse variance map to be used for inverse-variance filter

- *delensing.py*: delensing estimation pipeline. Out put will be all.pdf for CIB-lensing cross delensing and internal.pdf for internal delensing. The data for these curves in the pdf will be saved in folder ../primordial/.
- *product.py*: generate reconstruction klm and semi-analytical n0 nlkk.
- *klmplot.py*: delensing estimation pipeline. Out put will be all.pdf for CIB-lensing cross delensing and internal.pdf for internal delensing. The data for these curves in the pdf will be saved in folder ../primordial/.
