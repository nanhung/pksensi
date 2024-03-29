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
#' Pearce, R. G., Setzer, R. W., Strope, C. L., 
#' Wambaugh, J. F., & Sipes, N. S. (2017).
#' httk: R package for high-throughput toxicokinetics.
#' \emph{Journal of Statistical Software}, 79(4), 1-26.
#'
#' Zurlinden, T. J., & Reisfeld, B. (2016).
#' Physiologically based modeling of the pharmacokinetics of acetaminophen and
#' its major metabolites in humans using a Bayesian population approach.
#' \emph{European journal of drug metabolism and pharmacokinetics},
#'  41(3), 267-280.
#'
#' @examples
#' params <- c(F = 0.9, KA = 1.2, KE = 0.2, V = 1.5)
#' t <- seq(0, 12, 0.1)
#' C <-FFPK(params = params, time = t)
#' plot(t, C, type = "l", xlab = "time", ylab = "concentration")
#'
#' @export
FFPK <- function(params, time, dose = 1) { # nolint
  a <- (dose * params[1] * params[2]) / (params[4] * (params[2] - params[3]))
  c <- a * exp(- params[3] * time) - a * exp(- params[2] * time)
  return(c)
}

#' @export
#' @describeIn models Download pbtk1cpt.model file.
pbtk1cpt_model <- function() {
  mpath <- system.file("models", "pbtk1cpt.model", package = "pksensi")
  file.copy(mpath, getwd()) |> invisible()
}

#' @export
#' @describeIn models Download pbpk_apap.model file.
pbpk_apap_model <- function() {
  mpath <- system.file("models", "pbpk_apap.model", package = "pksensi")
  file.copy(mpath, getwd()) |> invisible()
}
