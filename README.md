# QSM Reconstruction Challenge 2016: Data and Code

Data and MATLAB scripts for the QSM Reconstruction Challenge 2016 (4th International Workshop on MRI Phase Contrast & QSM, Graz). The `data` folder contains magnitude, chi_33, COSMOS, brain mask, tissue phase and MPRAGE volumes; raw phase data are also included but are not used for the challenge.

`script_recon_public_v1.m` loads the single-orientation phase data and the chi_33 and COSMOS reference volumes, performs TKD and closed-form L2 dipole inversion, and evaluates the performance metrics. `utils` contains helper functions for display and evaluation. The challenge slides are included as `2016_05_05_Recon_Challenge_16by9.pdf`.

## Usage

Run `script_recon_public_v1.m` from the repository folder in MATLAB.

## Third-party code

`utils/compute_ssim.m` is based on the SSIM implementation by Zhou Wang (copyright notice retained in the file).

## Reference

C Langkammer et al. Quantitative susceptibility mapping: Report from the 2016 reconstruction challenge. *Magnetic Resonance in Medicine*, 2018. Challenge website: http://www.neuroimaging.at/pages/qsm.php

Contact: berkin AT nmr.mgh.harvard.edu
