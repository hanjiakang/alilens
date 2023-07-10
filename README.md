# alilens

Python 3.9.7

AliCPTlens is a python script for cosmology containning most of AliCPT CMB lensing pipeline, developed by Bin Hu, Jiakang Han, Jinyi Liu on behalf of the *AliCPT* collaboration ([publication here.](https://arxiv.org/abs/2303.05705)) This script is mostly based on [Plancklens](https://github.com/carronj/plancklens) by Julien Carron. 

This script may be used for:
* lensing reconstruction on *AliCPT* dataset consists of foregorund residue, CMB signal and instrumental noise
* delensing based on lensing reconstruction and its cross-correlation with Cosmic Microwave Background

### Requirements
To use alilens, you will need [Plancklens](https://github.com/carronj/plancklens), [lenspyx](https://github.com/carronj/lenspyx) and [healpy](https://github.com/carronj/lenspyx) installed locally. Check out their github for instructions. 

### Contents

This code contains most of the AliCPT lensing reconstruction pipeline (see Ali lensing paper [2022](https://arxiv.org/abs/2204.08158), [2023].(https://arxiv.org/abs/2303.05705)) Here we provide two versions corresponding to the above papers respectively:

* 2022 version consisting of data simulation scripts, which generate skymaps with planck FPP CMB signal and instrumental noise based on AliCPT noise variance map. It calculates the lensing reconstruciton results for both 95 GHz and 150 GHz from temperature only, polarization only and MV estimators.
* 2023 version consistiing of lensing reconstruction pipeline for polarization only estimation of pre-generated data base, which contains maps with CMB signal, instrumental noise and foreground residues after ILC fg removing method.

### 2022 version


### 2023 version
