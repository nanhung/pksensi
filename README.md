
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pksensi

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Dependencies](https://tinyverse.netlify.com/badge/pksensi)](https://cran.r-project.org/package=pksensi)
[![Travis-CI Build
Status](https://travis-ci.org/nanhung/pksensi.svg?branch=master)](https://travis-ci.org/nanhung/pksensi)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/nanhung/pksensi?branch=master&svg=true)](https://ci.appveyor.com/project/nanhung/pksensi)
[![Coverage
Status](https://codecov.io/gh/nanhung/pksensi/branch/master/graph/badge.svg)](https://codecov.io/gh/nanhung/pksensi?branch=master)  
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-last-release/pksensi)](https://cran.r-project.org/package=pksensi)
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/pksensi)](https://cran.r-project.org/package=pksensi)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/pksensi?color=orange)](https://CRAN.R-project.org/package=pksensi)
<!-- badges: end -->

**pksensi** implements the global sensitivity analysis workflow to
investigate the parameter uncertainty and sensitivity in physiologically
based kinetic (PK) models, especially the physiologically based
pharmacokinetic/toxicokinetic model with multivariate outputs. The
package also provides some functions to check the convergence and
sensitivity of model parameters.

Through **pksensi**, you can: - Run sensitivity analysis for PK models
in R with script that were written in C or GNU MCSim. - Decision
support: The output results and visualization tools can be used to
easily determine which parameters have “non-influential” effects on the
model output and can be fixed in model calibration.

## Installation

You can install the released version of pksensi from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pksensi")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("nanhung/pksensi")
```

-   This package includes a function that can help you install GNU MCsim
    more easily through the function `mcsim_install()`.

-   All updated details can be found in
    [NEWS.md](https://github.com/nanhung/pksensi/blob/master/NEWS.md).

-   **NOTE:** Windows users need to install
    [Rtools40](https://cran.r-project.org/bin/windows/Rtools/) to
    compile the model code.

## Workflow

![](https://i.ibb.co/tqpDLrk/sensitivity-workflow.png)

## Example

This is a basic example of applying **pksensi** in one-compartment pbtk
model:

``` r
library(pksensi)
## basic example code
```

# Step 1. Construct 1-cpt pbtk model

``` r
pbtk1cpt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dAgutlument = - kgutabs * Agutlument
    dAcompartment = kgutabs * Agutlument - ke * Acompartment
    dAmetabolized = ke * Acompartment
    Ccompartment = Acompartment / vdist * BW;
    list(c(dAgutlument, dAcompartment, dAmetabolized), 
         "Ccompartment" = Ccompartment) 
  })
}
```

# Step 2. Define initial conditions, output time steps and variable

``` r
initState <- c(Agutlument = 10, Acompartment = 0, Ametabolized = 0)
t <- seq(from = 0.01, to = 24.01, by = 1)
outputs <- c("Ccompartment")
```

# Step 3. Generate parameter matrix

## 3.1. (Optional) Extract parameter value from httk package

``` r
library(httk)
pars1comp <- (parameterize_1comp(chem.name = "acetaminophen"))
#> Human volume of distribution returned in units of L/kg BW.
```

## 3.2. Set parameter distributions

``` r
q <- c("qunif", "qunif", "qunif", "qnorm")
q.arg <- list(list(min = pars1comp$Vdist / 2, max = pars1comp$Vdist * 2),
              list(min = pars1comp$kelim / 2, max = pars1comp$kelim * 2),
              list(min = pars1comp$kgutabs / 2, max = pars1comp$kgutabs * 2),
              list(mean = pars1comp$BW, sd = 5))
```

## 3.3. Create parameter matrix

``` r
set.seed(1234)
params <- c("vdist", "ke", "kgutabs", "BW")
x <- rfast99(params, n = 200, q = q, q.arg = q.arg, replicate = 1)
```

# Step 4. Conduct simulation (will take few minutes)

``` r
out <- solve_fun(x, time = t, func = pbtk1cpt, initState = initState, outnames = outputs)
#> Starting time: 2021-06-04 15:50:47
#> Ending time: 2021-06-04 15:50:58
```

# Step 5. Uncertainty analysis

``` r
pksim(out)  # Use to compare with "real" data (if any)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

# Step 6. Check and visualize the result of sensitivity analysis

``` r
plot(out)   # Visualize result
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
check(out)  # Print result to console
#> 
#> Sensitivity check ( Index > 0.05 )
#> ----------------------------------
#> First order:
#>  vdist ke kgutabs 
#> 
#> Interaction:
#>  vdist ke kgutabs 
#> 
#> Total order:
#>  vdist ke kgutabs 
#> 
#> Unselected factors in total order:
#>  BW 
#> 
#> 
#> Convergence check ( Index > 0.05 )
#> ----------------------------------
#> First order:
#>   
#> 
#> Interaction:
#>   
#> 
#> Total order:
#> 
```

## Citation


    To cite pksensi in publications use:

      Hsieh, N-H., Reisfeld B., and Chiu W.A., (2018). pksensi: An R
      package to apply global sensitivity analysis in physiologically based
      kinetic modeling SoftwareX, 12, 100609.
      https://doi.org/10.1016/j.softx.2020.100609

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {{pksensi}: An R package to apply global sensitivity analysis in physiologically based kinetic modeling},
        author = {Nan-Hung Hsieh and Brad Reisfeld and Weihsueh A. Chiu},
        journal = {SoftwareX},
        year = {2020},
        volume = {12},
        pages = {100609},
        doi = {10.1016/j.softx.2020.100609},
      }

## Reference

Hsieh NH, Reisfeld B, Bois FY, Chiu WA. [Applying a global sensitivity
analysis workflow to improve the computational efficiencies in
physiologically-based pharmacokinetic
modeling](https://www.frontiersin.org/articles/10.3389/fphar.2018.00588/full).
Frontiers in Pharmacology 2018 Jun; 9:588.

Hsieh NH, Reisfeld B, Chiu WA. [pksensi: an R package to apply
sensitivity analysis in pharmacokinetic
modeling](https://nanhung.rbind.io/poster/2019-SOT.pdf). 58th SOT Annual
Meeting, Baltimore, USA, March 10–14, 2019.
