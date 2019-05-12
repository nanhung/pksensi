# pksensi 1.1.0

### Fix Bug:

* Used single time point in `solve_mcsim()`

### Update function:

**Installation**

* Updated `version = 6.1.0` in `install_mcsim()`
* Adopt to file name in MCSim_under_R - "model.R.exe"

**Modeling**

* Change assignment `n` to `monte_carlo` in `solve_mcsim()`
* Revise the default name of output to "sim.out" and "setpoint.out" in `solve_mcsim()`
* Added solving message to track time spend in `solve_mcsim()` 
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
