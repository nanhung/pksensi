#' Solved PK model with given parameters space
#'
#' @description
#' Download the latest or specific version of MCSim from the official website
#' \url{https://www.gnu.org/software/mcsim/} and install it to the system directory.
#'
#' @param version a character of MCSim version number.
#' @param directory a character to assign the directory to install the MCSim.
#' The default directory is under \code{/home/username} (Linux), \code{/Users/username} (MacOS), and C drive (Windows).
#'
#' @import getPass
#' @importFrom utils download.file
#'
#' @rdname install_mcsim
#' @export
install_mcsim = function(version = "6.0.1", directory = NULL) {
  if (.Platform$OS.type == "windows") {
    stop("The current function haven't supprot Windows system")
  } else {
    message("Start install")
    version<-version
    URL <- sprintf('http://ftp.gnu.org/gnu/mcsim/mcsim-%s.tar.gz', version)
    tf <- tempfile()
    download.file(URL, tf, mode = "wb")

    name <- Sys.info()[['user']]

    if (is.null(directory)){
      if (Sys.info()[['sysname']] == "Darwin"){
        exdir <- paste0("/Users/", name, sprintf('/mcsim-%s', version))
      } else if (Sys.info()[['sysname']] == "Linux") {
        exdir <- paste0("/home/", name, sprintf('/mcsim-%s', version))
      } else if (Sys.info()[['sysname']] == "Windows") {
        exdir <- paste0("C:", sprintf('/mcsim-%s', version))
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
    message(paste0("The MCSim " , sprintf('%s', version), " is installed under ", mcsim.directory))
    setwd(current.wd)
  }
}
