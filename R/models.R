#' Test PK Model for Sensitivity Analysis
#'
#' The example test model is Flip-flop kinetics (FFPK).
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
  A <- (dose * parameters[1])/( parameters[3] * ( parameters[1] - parameters[2]))
  C <- A * exp(- parameters[2] * time) - A * exp(- parameters[1] * time)
  return(C)
}
