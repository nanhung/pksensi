#' Extended Fourier Amplitude Sensitivity Test with Random Phase Shift
#'
#' @description
#'   \code{rfast99} is used to create the sequances for each parameter. It is based on the \code{fast99} function in \pkg{sensitivity} package.
#'
#' @importFrom stats runif fft var
#'
#' @param factors an integer for the giving number of factors, or a vector of character strings giving their names.
#' @param n an integer for the sampling number.
#' @param M an integer specifying the interference parameter. The default is 4.
#' @param omega a vector giving the set of frequencies.
#' @param q a vector of quantile functions names corresponding to wanted factors distributions.
#' @param q.arg a list of quantile functions parameters.
#' @param replicate an integer to define the number of replication. The default is 1.
#' @param conf the confidence level for replication confidence intervals. The default is 0.95.
#'
#' @importFrom graphics mtext
#' @rdname rfast99
#' @export
rfast99 <- function(factors, n, M = 4, omega = NULL,
                    q = NULL, q.arg = NULL, replicate = 1, conf = 0.95) {

  # factors numbers and names

  if (is.character(factors)) {
    X.labels <- factors
    p <- length(X.labels)
  } else {
    p <- factors
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
  } else if (FALSE %in% sapply(q.arg, is.list)) { # q.arg isn't a list of lists
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

  # object of class "fast18"

  x <- list(M = M, s = s, omega = omega, a = a, factors = factors,
            replicate = replicate, conf = conf, call = match.call())
  class(x) <- "rfast99"

  return(x)
}
