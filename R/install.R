#' Install MCSim
#'
#' Download the latest or specific version of GNU MCSim from the official website
#' (\url{https://www.gnu.org/software/mcsim/}) and install it to the system directory.
#'
#' This function aims to help users install GNU MCSim more easily.
#' However, if you can not install it through this function.
#' You might need to follow the instruction of GNU MCSim and install it,
#' manually: \url{https://www.gnu.org/software/mcsim/mcsim.html#Installation}
#'
#' The default \code{mxstp} is setting to 500.
#' The user can increase \code{mxstp} to avoid possible error return.
#' If you meet any error when conduct sensitivity analysis, you can this function to reinstall GNU MCSim and set the higher \code{mxstp}.
#' The default \code{directory} to install MCSim is under \code{/home/username} (Linux), \code{/Users/username} (MacOS),
#' and \code{C:/Users/} (windows). To install MCSim in Windows, be sure to install Rtools or MinGW first.
#'
#' @references
#' F.Y. Bois, and D. Maszle, 1997,
#' MCSim: A Monte Carlo Simulation Program,
#' \emph{Journal of Statistical Software}, 2(9): 1–60.
#'
#' @param version a character of MCSim version number.
#' @param directory a character to assign the directory to put the GNU MCSim files.
#' @param mxstep a numeric value to assign the maximum number of (internally defined) steps
#' allowed during one call to the solver.
#'
#' @import getPass
#' @importFrom utils download.file
#'
#' @examples
#' \dontrun{mcsim_install(version = 6.1.0, mxstep = 10000)}
#'
#' @references \url{https://www.gnu.org/software/mcsim/}
#'
#' @rdname mcsim
#'
#' @export
mcsim_install <- function(version = "6.1.0", directory = NULL, mxstep = 500) {

  #  if (.Platform$OS.type == "windows") {
  #    stop("The function haven't supprot Windows system")
  #  if (.Platform$OS.type == "windows") {
  #    if (!(devtools::find_rtools() == T)) {
  #      warning("The Rtools should be installed first")
  #    }}

  message("Start install")
  version<-version
  URL <- sprintf('http://ftp.gnu.org/gnu/mcsim/mcsim-%s.tar.gz', version)
  tf <- tempfile()
  download.file(URL, tf, mode = "wb")

  name <- Sys.info()[['user']]

  if (is.null(directory)){
    if (Sys.info()[['sysname']] == "Darwin"){
      exdir <- paste0("/Users/", name)
    } else if (Sys.info()[['sysname']] == "Linux") {
      exdir <- paste0("/home/", name)
    } else if (Sys.info()[['sysname']] == "Windows") {
      exdir <- paste0("c:/Users/", name)
    }
  } else {exdir <- directory}

  utils::untar(tf, exdir = exdir)

  current.wd <- getwd()

  if (is.null(directory)){
    if (Sys.info()[['sysname']] == "Darwin"){
      setwd(paste0("/Users/", name, sprintf('/mcsim-%s', version)))

      # The MacOS used clang as default compiler, the following command is used to switch to GCC
      Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep=";"))

    } else if (Sys.info()[['sysname']] == "Linux") {
      setwd(paste0("/home/", name, sprintf('/mcsim-%s', version)))
    } else if (Sys.info()[['sysname']] == "Windows") {
      setwd(paste0("c:/Users/", name, sprintf('/mcsim-%s', version)))
    }
  } else {setwd(paste0(directory, sprintf('/mcsim-%s', version)))}

  mcsim.directory <-getwd()

  if (mxstep != 500){
    file <- paste0(getwd(), "/sim/lsodes1.c")
    lsodes1.c <- readLines(file)
    new.mxstp0 <- paste0("mxstp0 = ", mxstep)
    mxstp0 <- gsub("mxstp0 = 500", new.mxstp0, lsodes1.c)
    cat(mxstp0, file=file, sep="\n")
  }

  if (.Platform$OS.type == "unix"){
    system("./configure")
    system("make")
    system("make check")

    input <- getPass::getPass("Authentication is required to install MCSim (Password): ")

    if (Sys.info()[['sysname']] == "Darwin"){
      system("sudo -kS make install", input=input)
    } else if (Sys.info()[['sysname']] == "Linux"){
      system("sudo -kS sh -c 'make install; ldconfig'", input=input)
    }
  } else if (.Platform$OS.type == "windows") {
    Sys.setenv(PATH = paste("c:/Rtools/mingw_64/bin", Sys.getenv("PATH"), sep=";"))
    Sys.setenv(PATH = paste("c:/MinGW/bin", Sys.getenv("PATH"), sep=";"))
    setwd(paste0(mcsim.directory,"/mod"))
    generate_config.h()
    setwd(paste0(mcsim.directory,"/sim"))
    generate_config.h()
    setwd(mcsim.directory)
    system(paste0("gcc -o ./mod/mod.exe ./mod/*.c"))
    if(file.exists("./mod/mod.exe")){
      cat(paste0("Created 'mod.exe'"))
    }
  }
  cat("\n")
  message(paste0("The MCSim " , sprintf('%s', version), " is installed. The sourced folder is under ", mcsim.directory))
  setwd(current.wd)
}


#' @export
#' @describeIn mcsim Return the version number of GNU MCSim.
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
    if (class(version) == "character"){
      message("The '", version, "' is found in ", exdir)
    }
  }
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
