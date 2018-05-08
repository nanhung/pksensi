# pksensi: Sensitivity Analysis for Pharmacokinetic Models

**pksensi** implements the global sensitivity analysis workflow to investigate the parameter sensitivity in pharmacokinetic (PK) models, especially the physiologically based pharmacokinetic (PBPK) model with multivariate output. The package also provide some functions to check the converge and sensitivity of model parameters.

This package can be installed via the devtools package using:  
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("nanhung/pksensi")
```
