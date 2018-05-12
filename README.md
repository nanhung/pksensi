# pksensi: Sensitivity analysis for pharmacokinetic model

[![Travis-CI Build Status](https://travis-ci.org/nanhung/pksensi.svg?branch=master)](https://travis-ci.org/nanhung/pksensi)
[![Coverage Status](https://img.shields.io/codecov/c/github/nanhung/pksensi/master.svg)](https://codecov.io/github/nanhung/pksensi?branch=master)


**pksensi** implements the global sensitivity analysis workflow to investigate the parameter uncertainty and sensitivity in pharmacokinetic (PK) models, especially the physiologically based pharmacokinetic (PBPK) model and advanced compartment absorption and transit (ACAT) model with multivariate output. The package also provide some functions to check the convergence and sensitivity of model parameters.

This package can be installed via the devtools package using:  
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("nanhung/pksensi")
```
