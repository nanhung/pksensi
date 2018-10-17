#' Example PK Model for Sensitivity Analysis
#'
#' The example test model is Flip-flop kinetics (FFPK).
#' \deqn{$C(t) = \frac{F \cdot D \cdot k_a}{(k_a-k_e)V} (e^{-k_et}-e^{-k_at})$/}
#'
#' @param time the given time-points.
#' @param parameters a parameter matrix containing the input sample.
#' @param dose a given dose.
#'
#' @name pk_model
#' @aliases pk_model
#' @rdname models
#' @export
FFPK <- function(parameters, time, dose = 1){
  A <- (dose * parameters[1] * parameters[2])/( parameters[4] * ( parameters[2] - parameters[3]))
  C <- A * exp(- parameters[3] * time) - A * exp(- parameters[2] * time)
  return(C)
}
