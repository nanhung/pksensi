#' Example PK Model for Sensitivity Analysis
#'
#' One of the example test model is flip-flop pharmacokinetic (FFPK) model. The time-dependent concentration can be written as:
#' \deqn{C(t) = \frac{F \cdot D \cdot k_a}{(k_a-k_e)V} (e^{-k_et}-e^{-k_at})}
#' where \eqn{F} is the fraction or percentage of the administrated dose that can reach the general circulation,
#' \eqn{k_a} is the first-order absorption rate constant (/time),
#' \eqn{k_e} is the first-order elimination rate constant (/time), and \eqn{V} is the distribution volume.
#'
#' @param time the given time-points.
#' @param params a parameter matrix containing the input sample.
#' @param dose a given dose.
#'
#' @name pk_model
#' @aliases pk_model
#' @rdname models
#'
#' @examples
#' params <- c(F = 0.9, KA = 1.2, KE = 0.2, V = 1.5)
#' t <- seq(0, 12, 0.1)
#' C <-FFPK(params = params, time = t)
#' plot(t, C, type = "l", xlab = "time", ylab = "concentration")
#'
#' @export
FFPK <- function(params, time, dose = 1){
  A <- (dose * params[1] * params[2])/( params[4] * ( params[2] - params[3]))
  C <- A * exp(- params[3] * time) - A * exp(- params[2] * time)
  return(C)
}

#' @export
#' @describeIn models Generate pbtk1comp.c file.
pbtk1comp.c = function(){
  url = "https://raw.githubusercontent.com/cran/httk/master/src/pbtk1comp.c"
  destfile = paste0(getwd(),"/pbtk1comp.c")
  download.file(url, destfile)
}

