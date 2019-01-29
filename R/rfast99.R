#' Extended Fourier Amplitude Sensitivity Test with Random Phase Shift
#'
#' \code{rfast99} is used to create the sequances for each parameter.
#' It is based on the \code{fast99} function in \pkg{sensitivity} package.
#'
#' @importFrom stats runif fft var
#'
#' @param params an integer for the giving number of parameters, or a vector of character strings giving their names.
#' @param n an integer for the sampling number.
#' @param M an integer specifying the interference parameter. The default is 4.
#' @param omega a vector giving the set of frequencies.
#' @param q a vector of quantile functions names corresponding to wanted parameters distributions.
#' @param q.arg a list of quantile functions parameters.
#' @param replicate an integer to define the number of replication. The default is set to 1 turn off the replication.
#' @param conf the confidence level for replication confidence intervals. The default is 0.95.
#'
#' @references
#' A. Saltelli, S. Tarantola and K. Chan, 1999, A quantitative, model independent method for
#' global sensitivity analysis of model output,
#' \emph{Technometrics}, 41, 39-56.
#'
#' R. I. Cukier, H. B. Levine and K. E. Schuler, 1978,
#' Nonlinear sensitivity analysis of multiparameter model systems.
#' \emph{Journal of Computational Physics}, 26, 1–42.
#'
#' @return
#' The returned parameter value will be stored in an array with
#' c(model evaluation, replication, parameters).
#'
#' @importFrom graphics mtext
#'
#' @examples
#' # Generate the parameter matrix with 20 replications
#' q <- "qunif"
#' q.arg <- list(min = 0, max = 1)
#'
#' set.seed(1234)
#' x <- rfast99(params = 3, n = 100, replicate = 20, q = q, q.arg = q.arg)
#' dim(x$a) # the array of c(model evaluation, replication, parameters).
#'
#' \dontrun{
#' save(x, file = "input_parameters.rda")
#' }
#'
#' @export
rfast99 <- function(params, n, M = 4, omega = NULL,
                    q = NULL, q.arg = NULL, replicate = 1, conf = 0.95) {

  # params numbers and names

  if (is.character(params)) {
    X.labels <- params
    p <- length(X.labels)
  } else {
    p <- params
    X.labels <- paste("X", 1 : p, sep = "")
  }

  # quantiles

  if (is.null(q)) {
    stop("Please assign the distribution(s) of quantile function")
  } else if (length(q) == 1) {
    q <- rep(q, p)
  }
  if (is.null(q.arg)) {
    stop("Please assign the arguments for defined distribution")
  } else if (FALSE %in% sapply(q.arg, is.list)) {
    q.arg <- rep(list(q.arg), p)
  }

  # set of frequencies

  if (is.null(omega)) {
    omega <- numeric(p)
    omega[1] <- floor((n - 1) / (2 * M))
    m <- floor(omega[1] / (2 * M))
    if (m >= p - 1) {
      omega[-1] <- floor(seq(from = 1, to = m, length.out = p - 1))
    } else {
      omega[-1] <- (0 : (p - 2)) %% m + 1
    }
  }

  # discretization of the s-space

  s <- 2 * pi / n * (0 : (n - 1))

  # transformation to get points in the x-space

  dim <- c(n * p, replicate, p)
  a <- array((n * p):( n * p * replicate * p), dim = dim)

  for (k in 1:replicate){
    X <- a[,k,]
    omega2 <- numeric(p)
    for (i in 1 : p) {
      omega2[i] <- omega[1]
      omega2[-i] <- omega[-1]
      l <- seq((i - 1) * n + 1, i * n)
      for (j in 1 : p) {
        v <- runif(1, min = 0, max = 2 * pi) # add random phase shift
        g <- 0.5 + 1 / pi * asin(sin(omega2[j] * s + v))
        X[l, j] <- do.call(q[j], c(list(p = g), q.arg[[j]]))
      }
    }
    # a[,k,] <- X
    a[,k,] <- round(X, digits = 4)
  }
  dimnames(a)[[3]] <- X.labels

  # object of class "rfast99"

  x <- list(M = M, s = s, omega = omega, a = a, params = params,
            replicate = replicate, conf = conf, call = match.call())
  class(x) <- "rfast99"

  return(x)
}
