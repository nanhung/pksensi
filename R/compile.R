#' Model Compiler
#'
#' The \code{compile_model} is used to compile the C file or MCSim's model file
#' to generate the executable file in numerical analysis.
#'
#' Generally, the solving function through MCSim can provide faster speed than exporting C in R.
#' Therefore, this function set \code{use_model_file = TRUE} and \code{application = 'mcsim'} as a default setting
#' and suggest to use MCSim to solve the differential equation.
#' To compile MCSim in Windows, be sure to install Rtools or MinGW first.
#' For Windows user, to compile MCSim's model file,
#' the \code{version} of MCSim should provide to conduct model compiling.
#'
#' @param mName a string giving the name of the model or C file (without extension).
#' @param use_model_file a logical value to operate the compiler to use model or C file,
#' the default is set to \code{TRUE} to assign the MCSim's model file in compiling.
#' @param application a character to assign the specific methods (\code{mcsim} or \code{R})
#' that will be applied to the numerical analysis (default is \code{mcsim}).
#' @param version a character to assign the version of MCSim that had been installed.
#' The version must be assigned for Windows user.
#'
#' @return
#' The default \code{application} is set to \code{'mcsim'}
#' to generate the executable file to solve differential equations by MCSim.
#' If \code{application = 'R'},
#' the function will compile and create dynamic-link libraries (.dll) on Windows and
#' shared objects (.so) on Unix-likes systems (e.g., Linux and MacOS).
#'
#'
#' @export
compile_model <- function (mName, application = 'mcsim', use_model_file = TRUE, version = NULL) {

  if (application == 'mcsim' && .Platform$OS.type == "windows"){
    mName <- paste0(mName,".model")
  }

#  if (.Platform$OS.type == "windows") {
#    if (!(devtools::find_rtools() == T)) {
#      warning("The Rtools should be installed first")
#    }
#  }

  if (use_model_file == T){ # Generate the ".c" file and "_inits.R" from model file
    if(file.exists(paste0(mName, ".model")) && .Platform$OS.type == "unix"){
      if (application == "mcsim"){
        system (paste0("mod ", mName, ".model ", mName,".c"))
      } else if (application == "R") {
        system (paste0("mod -R ", mName, ".model ", mName,".c")) # model to c file and include <R.h>
      } else {stop("Please assign the application to 'mcsim' or 'R'")}
    } else if (.Platform$OS.type == "windows") {
      if (is.null(version)) stop("Please provide the version of MCSim")
      mcsim <- paste0("mcsim-", version)
      name <- Sys.info()[['user']]
      exdir <- paste0("c:\\Users\\", name, "\\", mcsim)
      mod <-paste0("c:/Users/", name, "/", mcsim, "/mod/mod.exe")
      sim <-paste0("c:/Users/", name, "/", mcsim, "/sim")
    }
  }

  if (is.loaded("derivs", PACKAGE=mName)) dyn.unload(paste0(mName,.Platform$dynlib.ext))

  if (application == "mcsim"){
    if (.Platform$OS.type == "unix"){
      system(paste0("gcc -O3 -I/usr/local/include -L/usr/local/lib -g -O2 ", mName, ".c", " -lmcsim -o", " mcsim.", mName, " -lm -llapack -Wall"))

      exec <- paste0("mcsim.", mName)
      if (file.exists(exec)){
        cat(paste0("* Created executable file 'mcsim.", mName, "'."))
      } else stop("* Error in model compilation.")

    } else if ((.Platform$OS.type == "windows")) {
      Sys.setenv(PATH = paste("c:\\Rtools\\mingw_64\\bin", Sys.getenv("PATH"), sep=";"))
      Sys.setenv(PATH = paste("c:\\MinGW\\bin", Sys.getenv("PATH"), sep=";"))
      system(paste0(mod, " ", mName, " ", mName, ".c"))
      system(paste0("gcc -O3 -I.. -I", sim, " -o mcsim.", mName, ".exe ", mName, ".c ", sim, "/*.c", " -lm "))
      if (file.exists(paste0("mcsim.", mName, ".exe"))){
        cat(paste0("* Created executable file 'mcsim.", mName, ".exe'."))
      }
    }

  } else if (application == "R"){
    if (.Platform$OS.type == "windows" && use_model_file == T){
      system(paste0(mod, " -R ", mName, ".model ", mName, ".c"))
    }

    system (paste0("R CMD SHLIB ", mName, ".c")) # create .o and .so (or .dll) files
    if (file.exists(paste0(mName, ".so"))){
      cat(paste0("* Created file '", mName, ".so'."), "\n")
    }
    if (file.exists(paste0(mName, ".dll"))){
      cat(paste0("* Created file '", mName, ".dll'."), "\n")
    }
  }

  if (application == "R"){
    dyn.load(paste0(mName, .Platform$dynlib.ext))
  }

  if(file.exists(paste0(mName, "_inits.R"))){
    source(paste0(mName, "_inits.R"))
  }
}
