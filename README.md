
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pksensi <img src="man/figures/logo.png" height="200px" align="right" alt="pksensi logo" />

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/nanhung/pksensi/workflows/R-CMD-check/badge.svg)](https://github.com/nanhung/pksensi/actions)
[![Codecov test
coverage](https://codecov.io/gh/nanhung/pksensi/branch/master/graph/badge.svg)](https://codecov.io/gh/nanhung/pksensi?branch=master)
![CRAN downloads](https://cranlogs.r-pkg.org/badges/pksensi)
<!-- badges: end -->

**pksensi** implements the global sensitivity analysis workflow to
investigate the parameter uncertainty and sensitivity in physiologically
based kinetic (PK) models, especially the physiologically based
pharmacokinetic/toxicokinetic model with multivariate outputs. The
package also provides some functions to check the convergence and
sensitivity of model parameters.

Through **pksensi**, you can:

- Run sensitivity analysis for PK models in R with script that were
  written in C or GNU MCSim.

- Decision support: The output results and visualization tools can be
  used to easily determine which parameters have “non-influential”
  effects on the model output and can be fixed in following model
  calibration.

## Installation

You can install the released version of **pksensi** from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pksensi")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("nanhung/pksensi")
```

- This package includes a function that can help you install GNU MCSim
  more easily through the function `mcsim_install()`.

- All updated details can be found in
  [NEWS.md](https://github.com/nanhung/pksensi/blob/master/NEWS.md).

- **NOTE:** Windows users need to install
  [Rtools40](https://cran.r-project.org/bin/windows/Rtools/) to compile
  the model code.

## Workflow

<img src="man/figures/sensitivity-workflow.png" align="left" alt="Workflow of sensitivity analysis" />

**Note:** The parameter correlation (e.g., V<sub>max</sub> and
K<sub>M</sub> in metabolism) might be an issue in the global sensitivity
analysis. If you have experiment data, suggest using small datasets as a
sample in Markov Chain Monte Carlo Simulation. Then, check correlation
before conducting the sensitivity analysis. The issue will try to
address in the future version.

## Example

This is a basic example of applying **pksensi** in one-compartment pbtk
model:

``` r
library(pksensi)
```

### Step 1. Construct 1-cpt pbtk model

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

### Step 2. Define initial conditions, output time steps and variable

``` r
initState <- c(Agutlument = 10, Acompartment = 0, Ametabolized = 0)
t <- seq(from = 0.01, to = 24.01, by = 1)
outputs <- c("Ccompartment")
```

### Step 3. Generate parameter matrix

#### 3.1. (Optional) Extract parameter value from httk package

``` r
library(httk)
pars1comp <- (parameterize_1comp(chem.name = "acetaminophen"))
#> Human volume of distribution returned in units of L/kg BW.
```

#### 3.2. Set parameter distributions

``` r
q <- c("qunif", "qunif", "qunif", "qnorm")
q.arg <- list(list(min = pars1comp$Vdist / 2, max = pars1comp$Vdist * 2),
              list(min = pars1comp$kelim / 2, max = pars1comp$kelim * 2),
              list(min = pars1comp$kgutabs / 2, max = pars1comp$kgutabs * 2),
              list(mean = pars1comp$BW, sd = 5))
```

#### 3.3. Create parameter matrix

``` r
set.seed(1234)
params <- c("vdist", "ke", "kgutabs", "BW")
x <- rfast99(params, n = 200, q = q, q.arg = q.arg, replicate = 1)
```

### Step 4. Conduct simulation (will take few minutes with more replications)

``` r
out <- solve_fun(x, time = t, func = pbtk1cpt, initState = initState, outnames = outputs)
#> Starting time: 2024-11-27 07:07:16.85445
#> Ending time: 2024-11-27 07:07:25.711756
```

### Step 5. Uncertainty analysis

``` r
pksim(out)  # Use to compare with "real" data (if any)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" alt="pksim plot" width="100%" />

### Step 6. Check and visualize the result of sensitivity analysis

``` r
plot(out)   # Visualize result
```

<img src="man/figures/README-unnamed-chunk-9-1.png" alt="autoplot" width="100%" />

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

      Hsieh, N-H., Reisfeld B., and Chiu W.A., (2020). pksensi: An R
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

Hsieh NH, Reisfeld B, Chiu WA. [pksensi: An R package to apply global
sensitivity analysis in physiologically based kinetic
modeling](https://doi.org/10.1016/j.softx.2020.100609). SoftwareX 2020
Jul; 12:100609.

Hsieh NH, Bois FY, Tsakalozou E, Ni Z, Yoon M, Sun W, Klein M, Reisfeld
B, Chiu WA. [A Bayesian population physiologically based pharmacokinetic
absorption modeling approach to support generic drug development:
application to bupropion hydrochloride oral dosage
forms](https://doi.org/10.1007/s10928-021-09778-5). Journal of
Pharmacokinetics and Pharmacodynamics 2021 Sep; 22:1-6.
