# pksensi 1.1.7

* Suggested to use Rtools40
* Added parallel method in `solve_mcsim()`
* Fixed bug in `heat_check()`
* Added message of display info in `heat_check()` function
* Improved `heat_check()` to display the indices based on the setting cut-off

# pksensi 1.1.6

* Improved `heat_check()` to only show influential parameters as default
* Updated `check()` to return the SI and CI in a list

# pksensi 1.1.5

* `solve_mcsim()`: extense the flexible for the model in subfolder 
* Update the vignettes (add deSolve modeling example)

# pksensi 1.1.4

* Replace SE to NSE in `heat_check()`
* Update manual

# pksensi 1.1.3

* Update `solve_fun()` to run deSolve with Pure R code
* set default mxstep = 5000 in `mcsim_install()`

# pksensi 1.1.2

* Add timer in `solve_fun()`
* Add default MCSim `version = '6.1.0'` in `compile_model()`
* Update `pksim()`: Can customize xlab and ylab
* Update `pbtk1cpt model`: Add body weight
* Update vignettes
* Fix bug: Create `.dll` in win10 through `compile_model()`
* Fix bug: `plot()` cannot display variable name when using character
* Revise the `mcsim_version()` for Windows OS

# pksensi 1.1.1

### Fix Bug:

* Add `rtol` and `atol` in `solve_mcsim()`
* Switch default compiler from Clang to GCC in MacOS
* Fix the mName problem under Windows environment in `solve_mcsim()`

### Update vignette

* Remind to set the `mxstep` to 5000 in APAP-PBPK example
* Add note to APAP vignette (Set `atol` and `rtol`)

### Update function

* mcsim_version() can work in MCSim_under_R
* Stop the function when compile error in `compile_model()`
* Revise the default name of output to "simmc.out" and "setpts.out" in `solve_mcsim()`


# pksensi 1.1.0

### Fix Bug:

* Used single time point in `solve_mcsim()`

### Update function:

**Installation**

* Update the default MCSim version to `version = 6.1.0` in `install_mcsim()`
* Adopt the name of model file ("model.R.exe") to MCSim_under_R project

**Modeling**

* Change assignment `n` to `monte_carlo` in `solve_mcsim()`
* Revise the default name of output to "sim.out" and "setpoint.out" in `solve_mcsim()`
* Add solving message to track time spend in `solve_mcsim()` 
* Add assignment `tell = T` to automatically combine the output y in decoupling simulation x in `solve_fun()` and `solve_mcsim()`

**Plot**

* Transfer the log-transformed value to natural scale in `pksim()`
* Revise the discrete time condition to `length(times) < 16` in `heat_check()`


# pksensi 1.0.1

### New vignette:

* Add "PBTK 1-compartment model"
* Add "Acetaminophen-PBPK model"

### New example:

* Update example in `solve_mcsim()`

### New function:

* Add `pbtk1cpt_model()`
* Add `pbpk_apap_model()`
* Add function `mcsim_version()`

### NEW dataset:

* Add `APAP` dataset

### Fix Bug:

- Fix the `order` argument in `heat_check()`
- Add message in `generate_infile()`

### Update function:

* Remove argument `params` in `solve_fun()`
* No need to define `infile.name` and `outfile.name` in `solve_mcsim()`

### Change:

* Change function's name `install_mcsim()` to `mcsim_install()`


# pksensi 1.0.0

* Initial release in CRAN
