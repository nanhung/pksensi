#' Example PK Model for Sensitivity Analysis
#'
#' The example PK model that were used in sensitivity testing.
#' Three examples are included: Flip-flop pharmacokinetic model,
#' one-compartment toxicokinetic model from \pkg{httk} (Pearce et al. 2017),
#' and acetaminophen pharmacokinetic model (Zurlinden et al. 2016).
#'
#' @param time a numeric vector to define the given time point(s).
#' @param params a numeric vector to define the input parameter value.
#' @param dose a numeric value to define the given dose in flip-flop model.
#'
#' @name models
#' @aliases models
#' @rdname models
#'
#' @references
#' Pearce, R. G., Setzer, R. W., Strope, C. L., Wambaugh, J. F., & Sipes, N. S. (2017).
#' httk: R package for high-throughput toxicokinetics.
#' \emph{Journal of Statistical Software}, 79(4), 1-26.
#'
#' Zurlinden, T. J., & Reisfeld, B. (2016).
#' Physiologically based modeling of the pharmacokinetics of acetaminophen and
#' its major metabolites in humans using a Bayesian population approach.
#' \emph{European journal of drug metabolism and pharmacokinetics}, 41(3), 267-280.
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
#' @describeIn models Download pbtk1cpt.model file.
pbtk1cpt_model = function(){
  url = "https://raw.githubusercontent.com/nanhung/pksensi/master/tests/pbtk1cpt.model"
  destfile = paste0(getwd(),"/pbtk1cpt.model")
  download.file(url, destfile)
}

#' @export
#' @describeIn models Download pbpk_apap.model file.
pbpk_apap_model = function(){
  url = "https://raw.githubusercontent.com/nanhung/pksensi/master/tests/pbpk_apap.model"
  destfile = paste0(getwd(),"/pbpk_apap.model")
  download.file(url, destfile)
}
