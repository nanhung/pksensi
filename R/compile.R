#' Model Compiler
#'
#' The \code{compile_model} is used to compile the model code that is written under C or \pkg{GNU MCSim} format and
#' generate the executable program in numerical analysis.
#'
#' Generally, the solving function through \pkg{GNU MCSim} can provide faster computing speed than exporting C in R.
#' Therefore, this function set \code{use_model_file = TRUE} and \code{application = 'mcsim'}
#' as a default setting, suggesting to use \pkg{GNU MCSim} as main solver to solve the differential equation.
#' To compile model code in Windows, be sure to install Rtools (rtools40) first.
#' In addition, the \code{version} of \pkg{GNU MCSim} should provide to conduct model compiling in Windows.
#'
#' @param mName a string giving the name of the model code (without extension).
#' @param use_model_file a logical value to operate the compiler to model or C file,
#' the default is set to \code{TRUE} to assign the \pkg{GNU MCSim}'s model file in compiling.
#' @param application a character to assign the specific methods (\code{mcsim} or \code{R})
#' that will be applied to the numerical analysis (default is \code{mcsim}).
#' The version must be assigned for Windows user (default is \code{6.2.0}).
#' @param mcsim_dir a character to assign the location of MCSim directory
#' (The default is set to home).
#' @param mod_dir a character giving the name of the directory that is used to store the model file.
#'
#' @return
#' The default \code{application} is set to \code{'mcsim'}
#' to generate the executable program to solve differential equations by \pkg{GNU MCSim}.
#' If \code{application = 'R'},
#' the function will compile and create dynamic-link libraries (.dll) on Windows and
#' shared objects (.so) on Unix-likes systems (e.g., Linux and MacOS) that can link with \pkg{deSolve} solver.
#'
#' @export
compile_model <- function (mName, application = 'mcsim', use_model_file = TRUE,
                           mcsim_dir = NULL, mod_dir = NULL) {

  if (application == 'mcsim' && .Platform$OS.type == "windows"){
    mName <- paste0(mName,".model")
  }

  if (.Platform$OS.type == "unix"){
    if (is.null(mcsim_dir)) mcsim_directory <- paste0(Sys.getenv("HOME"), "/mcsim")
    else mcsim_directory <- paste0(mcsim_dir, "/mcsim")

    bin_path <- paste0(mcsim_directory, "/bin")
    Sys.setenv(PATH = paste(bin_path, Sys.getenv("PATH"), sep=":"))
    Sys.setenv(LD_LIBRARY_PATH = paste(mcsim_directory, "/lib",
                                       Sys.getenv("LD_LIBRARY_PATH"), sep=":"))
  }

  if (use_model_file == T){ # Generate the ".c" file and "_inits.R" from model file
    if(file.exists(paste0(mName, ".model")) && .Platform$OS.type == "unix"){
      if (application == "mcsim"){
        system (paste0("mod ", mName, ".model ", mName,".c"))
      } else if (application == "R") {
        system (paste0("mod -R ", mName, ".model ", mName,".c")) # model to c file and include <R.h>
      } else {stop("Please assign the application to 'mcsim' or 'R'")}
    } else if (.Platform$OS.type == "windows") {

      home_dir <- Sys.getenv("HOME")
      if(is.null(mcsim_dir)) mcsim_dir <- paste0(home_dir, "/mcsim") else
        mcsim_dir <- paste0(mcsim_dir, "/mcsim")
      mod <- paste0(mcsim_dir, "/mod.exe")
      sim <-paste0(mcsim_dir, "/sim")
    }
  }

  if (is.loaded("derivs", PACKAGE=mName)) dyn.unload(paste0(mName,.Platform$dynlib.ext))

  if (application == "mcsim"){
    if (.Platform$OS.type == "unix"){
      #system(paste0("gcc -O3 -I/usr/local/include -L/usr/local/lib -g -O2 ", mName, ".c", " -lmcsim -o", " mcsim.", mName, " -lm -llapack -Wall"))

      if (is.null(mod_dir)) makemcsim <- paste0("makemcsims ", mName, ".model") else
        makemcsim <- paste0("makemcsims ", mod_dir, "/", mName, ".model")
      system(makemcsim)

      exec <- paste0("mcsim.", mName)
      if (!file.exists(exec)) stop("* Error in model compilation.\n")

    } else if ((.Platform$OS.type == "windows")) {
      Sys.setenv(PATH = paste("C:\\rtools40\\mingw64\\bin", Sys.getenv("PATH"), sep=";"))
      Sys.setenv(PATH = paste("C:\\rtools40\\bin", Sys.getenv("PATH"), sep=";"))

      if (is.null(mod_dir)) mod_c <- paste0(mod, " ", mName, " ", mName, ".c") else
        mod_c <- paste0(mod, " ", mod_dir, "/", mName, " ", mName, ".c")
      message("Creating model.c file...")
      system(mod_c)
      message("Creating modeling program...")
      makemcsim <- paste0("gcc -O3 -I.. -I", sim, " -o mcsim.", mName, ".exe ", mName, ".c ", sim, "/*.c", " -lm ")
      system(makemcsim)

      if (file.exists(paste0("mcsim.", mName, ".exe"))) message("Done.") else
        message ("Compiling error.")
    }

  } else if (application == "R"){
    if (.Platform$OS.type == "windows" && use_model_file == T){
      system(paste0(mod, " -R ", mName, ".model ", mName, ".c"))
      Sys.setenv(PATH = paste("C:\\rtools40\\mingw64\\bin", Sys.getenv("PATH"), sep=";"))
      Sys.setenv(PATH = paste("C:\\rtools40\\bin", Sys.getenv("PATH"), sep=";"))
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
