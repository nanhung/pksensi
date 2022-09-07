#' Download and install \pkg{GNU MCSim}
#'
#' Download the latest or specific version of \pkg{GNU MCSim} from the official website
#' (\url{https://www.gnu.org/software/mcsim/}) and install it to the system directory.
#'
#' This function aims to help users download (source: https://ftp.gnu.org/gnu/mcsim/)
#' and install \pkg{GNU MCSim} more easily.
#' However, if you can not install it through this function.
#' The additional way is to follow the instruction and install it manually:
#' \url{https://www.gnu.org/software/mcsim/mcsim.html#Installation}
#'
#' The default \code{mxstp} is setting to 5000.
#' The user can increase \code{mxstp} to avoid possible error return.
#' If you meet any error when conduct sensitivity analysis,
#' you can use this function to re-install \pkg{GNU MCSim} and set the higher \code{mxstp}.
#' The default installed \code{directory} is under \code{/home/username} (Linux),
#' \code{/Users/username} (MacOS),
#' and \code{C:/Users/username/Documents} (windows).
#' To install \pkg{GNU MCSim} in Windows, be sure to install Rtools first.
#' The current suggested Rtools version is 4.0.
#'
#' @references
#' Bois, F. Y., & Maszle, D. R. (1997).
#' MCSim: a Monte Carlo simulation program.
#' \emph{Journal of Statistical Software}, 2(9): 1–60.
#'
#' @param version a character of version number.
#' @param install_dir a character to assign the installed directory.
#' @param mxstep a numeric value to assign the maximum number of (internally defined) steps
#' allowed during one call to the solver.
#'
#' @importFrom utils download.file
#'
#' @references \url{https://www.gnu.org/software/mcsim/}
#'
#' @rdname mcsim
#'
#' @export
mcsim_install <- function(version = "6.2.0", install_dir = NULL, mxstep = 5000) {

  if (Sys.info()[['sysname']] == "Windows") {
    if(Sys.which("gcc") == "") stop("Please check the installation of Rtools.")
  } else if (Sys.info()[['sysname']] == "Linux") {
    if(Sys.which("gcc") == "") stop("Please check the installation of gcc.")
  }

  message("Start install")
  version<-version
  URL <- sprintf('http://ftp.gnu.org/gnu/mcsim/mcsim-%s.tar.gz', version)
  tf <- tempfile()
  download.file(URL, tf, mode = 'wb')

  name <- Sys.info()[['user']]
  home_dir <- Sys.getenv("HOME")

  # Defined directory (exdir) to place mcsim source code
  if (is.null(install_dir)){
    if (Sys.info()[['sysname']] == "Darwin"){
      exdir <- paste0("/Users/", name)
    } else if (Sys.info()[['sysname']] == "Linux") {
      exdir <- paste0("/home/", name)
    } else if (Sys.info()[['sysname']] == "Windows") {
      exdir <- home_dir
    }
  } else {exdir <- install_dir}

  utils::untar(tf, exdir = exdir)

  current.wd <- getwd() # the current working directory

  # Defined MCSim directory
  if (is.null(install_dir)){
    if (Sys.info()[['sysname']] == "Darwin"){
      setwd(paste0("/Users/", name, sprintf('/mcsim-%s', version)))

      # The MacOS used clang as default compiler, the following command is used to switch to GCC
      Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep=":"))

    } else if (Sys.info()[['sysname']] == "Linux") {
      setwd(paste0("/home/", name, sprintf('/mcsim-%s', version)))
    } else if (Sys.info()[['sysname']] == "Windows") {
      setwd(paste0(home_dir, sprintf('/mcsim-%s', version)))
    }
  } else {setwd(paste0(install_dir, sprintf('/mcsim-%s', version)))}

  mcsim.directory <- getwd()

  if (mxstep != 500){
    file <- paste0(getwd(), "/sim/lsodes1.c")
    lsodes1.c <- readLines(file)
    new.mxstp0 <- paste0("mxstp0 = ", mxstep)
    mxstp0 <- gsub("mxstp0 = 500", new.mxstp0, lsodes1.c)
    cat(mxstp0, file=file, sep="\n")
  }

  # Defined mcsim directory to place compiled files (e.g., bin, lib, etc)
  if (.Platform$OS.type == "unix"){

    if (is.null(install_dir)) mcsim_dir <- paste0(home_dir, "/mcsim") else
      mcsim_dir <- paste0(install_dir, "/mcsim")
    if(!dir.exists(mcsim_dir)) dir.create(mcsim_dir)

    system(paste0("./configure prefix=", mcsim_dir))
    system("make")
    system("make install")
    system("make check")

    bin_path <- paste0(mcsim_dir, "/bin")
    Sys.setenv(PATH = paste(bin_path, Sys.getenv("PATH"), sep=":"))
    lib_path <- paste0(mcsim_dir, "/lib")
    Sys.setenv(LD_LIBRARY_PATH = paste(lib_path,
                                       Sys.getenv("LD_LIBRARY_PATH"), sep=":"))

    message("\nChecking...")
    cat(Sys.which("makemcsim"))

    makemcsim <- paste0(bin_path, "/makemcsim")
    if(Sys.which("makemcsim") == makemcsim)
      message(paste0("\nThe MCSim " , sprintf('%s', version), " is installed."))
      message(paste0("The sourced folder is under ", mcsim.directory))

  } else if (.Platform$OS.type == "windows") {

    if(Sys.which("gcc") == ""){ # echo $PATH
      # Suggest use Rtools40
      PATH = "C:\\rtools40\\mingw64\\bin; C:\\rtools40\\usr\\bin"
      Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))

    } # PATH=$PATH:/c/Rtools/mingw_32/bin; export PATH

    setwd(paste0(mcsim.directory,"/mod"))
    generate_config.h()

    mcsim_sim <- paste0(mcsim.directory,"/sim")
    setwd(mcsim_sim)
    generate_config.h()
    sim_files <- list.files(mcsim_sim)

    # Create 'mcsim' directory
    if (is.null(install_dir)) {
      mcsim_dir <- paste0(home_dir, "/mcsim")
    } else mcsim_dir <- paste0(install_dir, "/mcsim")

    # Remove the previous installed version
    if(dir.exists(mcsim_dir)) {
      system(paste0("rm -rf ", mcsim_dir))
    #  if (menu(c("Yes", "No"),
    #           title = paste0("\nThe 'mcsim' directory is existed. ",
    #                          "Do you want to replace it?")) == 1)
    #    system(paste0("rm -rf ", mcsim_dir))
    #  else return(invisible())
    }

    dir.create(mcsim_dir)

    mcsim_sim_dir <- paste0(mcsim_dir, "/sim")
    if(!dir.exists(mcsim_sim_dir)) dir.create(mcsim_sim_dir)


    file.copy(from = paste0(mcsim_sim, "/", sim_files),
              to = paste0(mcsim_sim_dir, "/", sim_files))
    setwd(mcsim.directory)

    mod <- paste0(mcsim_dir, "/mod.exe")
    system(paste0("gcc -o ",  mod, " ./mod/*.c"))
    Sys.setenv(PATH = paste(mcsim_dir, Sys.getenv("PATH"), sep=";"))

    message("\nChecking...")
    cat(Sys.which("mod"))

    if(file.exists(mod)){
      message(paste0("\nThe MCSim " , sprintf('%s', version), " is installed."))
      message(paste0("The sourced folder is under ", mcsim.directory))
    } else stop("Cannot find mod.exe.")

  }
  cat("\n")
  setwd(current.wd)
}

#' @export
#' @describeIn mcsim Return the version number.
mcsim_version <- function(){

  if (.Platform$OS.type == "unix"){
    if(file.exists("MCSim/mod.exe")){ # for MCSim under R
      invisible(system("./MCSim/mod.exe -h | tee mod.mcsim.txt", intern = TRUE))
    } else {
      invisible(system("mod -h | tee mod.mcsim.txt", intern = TRUE))
    }
    l <- readLines("mod.mcsim.txt")
    invisible(file.remove("mod.mcsim.txt"))
    version <- substr(l[4], 6, 10)
    message("The current GNU MCSim version is ", version)
  }

  if (.Platform$OS.type == "windows") {
    name <- Sys.info()[['user']]
    exdir <- paste0("c:/Users/", name)
    l <- list.files(path = exdir)
    version <- l[grep("mcsim", l)]
    if (is.character(version)){
      message("The '", version, "' is found in ", exdir)
    }
  }
}

#' @export
#' @describeIn mcsim the function to set Rtools40 path
set_rtools40_path <- function(){

  if (.Platform$OS.type == "windows"){
    PATH = "C:\\rtools40\\mingw64\\bin; C:\\rtools40\\usr\\bin"
    Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))

    env.file <- ".Renviron"
    if (!file.exists(env.file)) file.create(env.file)
    path <- 'PATH="${RTOOLS40_HOME}\\usr\\bin;${RTOOLS40_HOME}\\mingw64\\bin;${PATH}"'
    write(path, file = paste0(getwd(), "/.Renviron"), append = TRUE)
  } else stop("You are not using Windows OS")

}

generate_config.h <- function(){
  cat("#define HAVE_DLFCN_H 1 \n",
      "#define HAVE_ERFC 1 \n",
      "#define HAVE_FLOAT_H 1 \n",
      "#define HAVE_FLOOR 1 \n",
      "#define HAVE_INTTYPES_H 1 \n",
      #"#define HAVE_LIBGSL 1 \n",
      "#define HAVE_LIBGSLCBLAS 1 \n",
      "#define HAVE_LIBLAPACK 1 \n",
      "#define HAVE_LIBM 1 \n",
      #"#define HAVE_LIBSBML 1 \n",
      #"#define HAVE_LIBSUNDIALS_CVODES 1 \n",
      #"#define HAVE_LIBSUNDIALS_NVECSERIAL 1 \n",
      "#define HAVE_LIMITS_H 1 \n",
      "#define HAVE_MALLOC 1 \n",
      "#define HAVE_MEMORY_H 1 \n",
      "#define HAVE_MODF 1 \n",
      "#define HAVE_POW 1 \n",
      "#define HAVE_REALLOC 1 \n",
      "#define HAVE_SQRT 1 \n",
      "#define HAVE_STDINT_H 1 \n",
      "#define HAVE_STDLIB_H 1 \n",
      "#define HAVE_STRCHR 1 \n",
      "#define HAVE_STRINGS_H 1 \n",
      "#define HAVE_STRING_H 1 \n",
      "#define HAVE_SYS_STAT_H 1 \n",
      "#define HAVE_SYS_TYPES_H 1 \n",
      "#define HAVE_UNISTD_H 1 \n",
      file = "config.h",
      sep = "")
}

set_permanent_path <- function(mcsim_dir=NULL){

  if (.Platform$OS.type == "unix"){

    if (Sys.getenv("SHELL")=="/usr/bin/zsh"){
      zshrc <- paste0(Sys.getenv("HOME"), "/.zshrc")
      cat('\n# Path to MCSim\n',
          file = zshrc, append=TRUE)
      if (is.null(mcsim_dir)){
        cat('export PATH=$PATH:$HOME/mcsim/bin',
            file = zshrc, append=TRUE)
      } else cat(paste0('export PATH=$PATH:', mcsim_dir, '/mcsim/bin'),
                 file = zshrc, append=TRUE)
      message(paste0("Done. Please check the .zshrc file at ", zshrc))
      } else if (Sys.getenv("SHELL")=="/bin/bash") {
        bashrc <- paste0(Sys.getenv("HOME"), "/.bashrc")
        cat('\n# Path to MCSim\n',
            file = bashrc, append=TRUE)
        if (is.null(mcsim_dir)){
          cat('export PATH=$PATH:$HOME/mcsim/bin',
              file = bashrc, append=TRUE)
        } else cat(paste0('export PATH=$PATH:', mcsim_dir, '/mcsim/bin'),
                   file = bashrc, append=TRUE)
        message(paste0("Done. Please check the .bashrc file at ", zshrc))
        }
  } else stop("The function is design for 'unix' system")
}
