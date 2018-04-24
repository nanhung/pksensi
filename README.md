# rfast99: Extended Fourier Amplitude Sensitivity Test with Random Phase Shift

**rfast99** implements the extended Fourier Amplitude Sensitivity Test (eFAST) method [(Saltelli et al. 1999)](https://www.tandfonline.com/doi/abs/10.1080/00401706.1999.10485594) with random phase shift to investigate the converge of sensitivity indices of model parameters under the given sampling number. The base codes were sourced from R [**sensitivity**](https://cran.r-project.org/web/packages/sensitivity/index.html) package.

This package can be installed via the devtools package using:  
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("nanhung/rfast99")
```
