#' Install MCSim
#'
#' @description
#' Download the latest or specific version of MCSim from the official website
#' \url{https://www.gnu.org/software/mcsim/} and install it to the system directory.
#' Currently, this function only support Unix-based systems (MacOS and Linux).
#'
#' @param version a character of MCSim version number.
#' @param directory a character to assign the directory to put the MCSim sourced folder.
#' The default directory is under \code{/home/username} (Linux), \code{/Users/username} (MacOS), and C drive (Windows).
#' @param mxstep a numeric value to assign the maximum number of (internally defined) steps
#' allowed during one call to the solver (default is 500). The user may increase mxstep to avoid this error return.
#'
#' @import getPass
#' @importFrom utils download.file
#'
#' @rdname install_mcsim
#' @export
install_mcsim = function(version = "6.0.1", directory = NULL, mxstep = 500) {
  if (.Platform$OS.type == "windows") {
    stop("The function haven't supprot Windows system")
  } else {
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
        exdir <- paste0("C:")
      }
    } else {exdir <- directory}

    utils::untar(tf, exdir = exdir)

    current.wd <- getwd()

    if (is.null(directory)){
      if (Sys.info()[['sysname']] == "Darwin"){
        setwd(paste0("/Users/", name, sprintf('/mcsim-%s', version)))
      } else if (Sys.info()[['sysname']] == "Linux") {
        setwd(paste0("/home/", name, sprintf('/mcsim-%s', version)))
      } else if (Sys.info()[['sysname']] == "Windows") {
        setwd(paste0("C:/", name, sprintf('/mcsim-%s', version)))
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

    system("./configure")
    system("make")
    system("make check")

    input <- getPass::getPass("Authentication is required to install MCSim (Password): ")

    if (Sys.info()[['sysname']] == "Darwin"){
      system("sudo -kS make install", input=input)
    } else if (Sys.info()[['sysname']] == "Linux"){
      system("sudo -kS sh -c 'make install; ldconfig'", input=input)
    }

    cat("\n")
    message(paste0("The MCSim " , sprintf('%s', version), " is installed. The sourced folder is under ", mcsim.directory))
    setwd(current.wd)
  }
}
