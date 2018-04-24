# rfast99: Extended Fourier Amplitude Sensitivity Tes with random randon phase shift

**rfast99** implements the "extended-FAST" method (Saltelli et al. 1999) with random phase shift to investigate the converge of sensitivity indices of model parameters under the given sampling number. The base codes were sourced from R [**sensitivity**](https://cran.r-project.org/web/packages/sensitivity/index.html) package.

This package can be installed via the devtools package using:  
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("nanhung/rfast99")
```
