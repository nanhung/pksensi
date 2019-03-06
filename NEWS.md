# pksensi 1.1.0

### Fix Bug:

* Used single time point in `solve_mcsim()`

### Update function:

* Added solving message in `solve_mcsim()`
* Updated `version = 6.1.0` in `install_mcsim()`

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
