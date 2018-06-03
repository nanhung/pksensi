#' Model Compiler
#'
#' @description
#' This function is used to compile the model file or C file to generate the executable file in numerical analysis
#'
#' @param mName a string giving the name of the model or C file (without extension).
#' @param model a logical value to operate the compiler to use model or C file,
#' the default is set to FALSE to assign the C file in compiling.
#' @param application a character to assign the specific methods (\code{mcsim} or \code{R})
#' that will be applied to the numerical analysis (default is \code{mcsim}).
#'
#' @importFrom devtools find_rtools
#'
#' @rdname compile
#' @export
compile <- function (mName, model = F, application = 'mcsim') {
  if (model == T){ # Generate the ".c" file and "_inits.R" from model file
    if(file.exists(paste0(mName, ".model")) && .Platform$OS.type == "unix"){
      if (application == "mcsim"){
        system (paste0("mod ", mName, ".model ", mName,".c"))
      } else if (application == "R") {
        system (paste0("mod -R ", mName, ".model ", mName,".c")) # model to c file and include <R.h>
      } else {stop("Please assign the application to 'mcsim' or 'R'")}
    } else if (.Platform$OS.type == "windows") {
      stop("The current function can't generate c file in Windows system")
    }
  }

  if (is.loaded("derivs", PACKAGE=mName)) dyn.unload(paste0(mName,.Platform$dynlib.ext))

  if (.Platform$OS.type == "windows") {
    if (!(devtools::find_rtools() == T)) {
      warning("The Rtools should be installed first")
    }
  }

  if (application == "mcsim"){
    system(paste0("gcc -O3 -I/usr/local/include -L/usr/local/lib -g -O2 ", mName, ".c", " -lmcsim -o", " mcsim.", mName, " -lm -llapack -Wall"))
    cat(paste0("* Created executable file 'mcsim.", mName, "'."))
  } else if (application == "R"){
    system (paste0("R CMD SHLIB ", mName, ".c")) # create .o and .so (or .dll) files
    cat(paste0("\n* Created executable file '", mName, ".so' or '", mName, ".dll'."))
  }


  if (application == "R"){
    dyn.load(paste0(mName, .Platform$dynlib.ext))
  }

  if(file.exists(paste0(mName, "_inits.R"))){
    source(paste0(mName, "_inits.R"))
  }
}
