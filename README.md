# pksensi <img src="man/figures/logo.png" height="200px" align="right" />

[![Travis-CI Build Status](https://travis-ci.org/nanhung/pksensi.svg?branch=master)](https://travis-ci.org/nanhung/pksensi)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/nanhung/pksensi?branch=master&svg=true)](https://ci.appveyor.com/project/nanhung/pksensi)
[![Coverage Status](https://img.shields.io/codecov/c/github/nanhung/pksensi/master.svg)](https://codecov.io/github/nanhung/pksensi?branch=master)  
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-last-release/pksensi)](https://cran.r-project.org/package=pksensi)
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/pksensi)](https://cran.r-project.org/package=pksensi)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/pksensi?color=orange)](https://CRAN.R-project.org/package=pksensi)

**pksensi** implements the global sensitivity analysis workflow to investigate the parameter uncertainty and sensitivity in pharmacokinetic (PK) models, especially the physiologically based pharmacokinetic (PBPK) model with multivariate outputs. The package also provide some functions to check the convergence and sensitivity of model parameters.

```
# To get pksensi from CRAN (current version 1.0.0):
install.packages("pksensi")

# Or get the the development version from GitHub:
if (!require("devtools")) install.packages("devtools")
devtools::install_github("nanhung/pksensi")
```

## Reference

Hsieh NH, Reisfeld B, Bois FY, Chiu WA. [Applying a global sensitivity analysis workflow to improve the computational efficiencies in physiologically-based pharmacokinetic modeling](https://www.frontiersin.org/articles/10.3389/fphar.2018.00588/full). Frontiers in Pharmacology 2018 Jun; 9:588.
