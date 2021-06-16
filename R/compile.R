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
#' @param version a character to assign the version of \pkg{GNU MCSim} that had been installed.
#' The version must be assigned for Windows user (default is \code{6.2.0}).
#'
#' @return
#' The default \code{application} is set to \code{'mcsim'}
#' to generate the executable program to solve differential equations by \pkg{GNU MCSim}.
#' If \code{application = 'R'},
#' the function will compile and create dynamic-link libraries (.dll) on Windows and
#' shared objects (.so) on Unix-likes systems (e.g., Linux and MacOS) that can link with \pkg{deSolve} solver.
#'
#' @export
compile_model <- function (mName, application = 'mcsim', use_model_file = TRUE, version = '6.2.0') {

  if (application == 'mcsim' && .Platform$OS.type == "windows"){
    mName <- paste0(mName,".model")
  }

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
      exdir <- paste0("c:/Users/", name, "/", mcsim)
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
        cat(paste0("* Created executable file 'mcsim.", mName, "'.\n"))
      } else stop("* Error in model compilation.\n")

    } else if ((.Platform$OS.type == "windows")) {
      Sys.setenv(PATH = paste("C:\\rtools40\\mingw64\\bin", Sys.getenv("PATH"), sep=";"))
      Sys.setenv(PATH = paste("C:\\rtools40\\bin", Sys.getenv("PATH"), sep=";"))
      system(paste0(mod, " ", mName, " ", mName, ".c"))
      system(paste0("gcc -O3 -I.. -I", sim, " -o mcsim.", mName, ".exe ", mName, ".c ", sim, "/*.c", " -lm "))
      if (file.exists(paste0("mcsim.", mName, ".exe"))){
        cat(paste0("* Created executable file 'mcsim.", mName, ".exe'."))
      }
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

compile_model_pkg <- function(mName, application = 'mcsim', version = '6.2.0'){

  if (.Platform$OS.type == "unix"){
    exe_file <- paste0("mcsim.", mName)
  } else if (.Platform$OS.type == "windows") exe_file <- paste0("mcsim.", mName, ".exe")

  mcsimdir <- system.file("mcsim", package = "pksensi")
  moddir <- paste0(mcsimdir, "/mcsim-", version, "/mod")
  simdir <- paste0(mcsimdir, "/mcsim-", version, "/sim")
  modpath <- paste0(moddir, "/mod.exe")

  if (file.exists(modpath) == F) mcsim_pkg(version = version)

  if (application == 'R'){
    system(paste0(modpath, " -R ", mName, ".model ", mName, ".c"))
    system (paste0("R CMD SHLIB ", mName, ".c")) # create *.dll files
    dyn.load(paste(mName, .Platform$dynlib.ext, sep="")) # load *.dll
    source(paste0(mName,"_inits.R"))
  } else if (application == 'mcsim') {
    system(paste0(modpath, " ", mName, ".model ", mName, ".c"))
    message(paste0("Compiling..."))
    system(paste0("gcc -O3 -I.. -I", simdir, " -o ", exe_file, " ", mName, ".c ", simdir, "/*.c -lm "))
    invisible(file.remove(paste0(mName, ".c")))
    if(file.exists(exe_file)) message(paste0("* Created executable program '", exe_file, "'."))
  }
}


