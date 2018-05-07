# Sobol function
sobol.fun<-function (parameters) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * parameters[j] - 2) + a[j])/(1 + a[j])
  }
  y
}

# Flip-Flop Kinetics
FFPK <- function(parameters, times, dose = 1){
  A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
  C <- A * exp(- parameters[2] * times) - A * exp(- parameters[1] * times)
  return(C)
}
