# Flip-Flop Kinetics
FFPK <- function(parameters, times, dose = 1){
  A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
  C <- A * exp(- parameters[2] * times) - A * exp(- parameters[1] * times)
  return(C)
}
