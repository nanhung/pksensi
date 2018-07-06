#' Extended Fourier Amplitude Sensitivity Test with Random Phase Shift
#'
#' @description
#'   \code{rfast99} is based on the fast99 function in \pkg{sensitivity} package.
#'
#' @importFrom stats runif fft var
#'
#' @param factors an integer for the giving number of factors, or a vector of character strings giving their names.
#' @param n an integer for the sampling number.
#' @param M an integer specifying the interference parameter.
#' @param omega a vector giving the set of frequencies.
#' @param q a vector of quantile functions names corresponding to wanted factors distributions.
#' @param q.arg a list of quantile functions parameters.
#' @param replicate an integer to define the number of replication, the default is 1.
#' @param conf the confidence level for replication confidence intervals.
#'
#' @details
#' The \code{print} method can used to print sensitivity and convergence indices
#' with given time-step. The \code{check} method provide the summary of
#' parameter sensitivity and convergence. The \code{plot} function provide the
#' time-course functional output of both indices for each parameter.
#' @param x a list of class "rfast99"
#' @param digits a integer to rounds the values in its first argument.
#' to the specified number of decimal places (default 4).
#' @param SI a numeric value to set the cut-off point for sensitivity index (default 0.1).
#' @param CI a numeric vlaue to set the cut-off point for convergence index (default 0.01).
#' @param ... additional arguments to customize the graphical parameters.
#'
#' @importFrom graphics mtext
#' @rdname rfast99
#' @export
rfast99 <- function(factors, n, M = 4, omega = NULL,
                    q = NULL, q.arg = NULL, replicate = 1, conf = 0.95, ...) {

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
    q <- rep("qunif", p)
  } else if (length(q) == 1) {
    q <- rep(q, p)
  }
  if (is.null(q.arg)) {
    q.arg <- rep(list(), p)
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
    a[,k,] <- X
  }
  dimnames(a)[[3]] <- X.labels

  # object of class "fast18"

  x <- list(M = M, s = s, omega = omega, a = a, factors = factors,
            replicate = replicate, conf = conf, call = match.call())
  class(x) <- "rfast99"

  return(x)
}

#' @method print rfast99
#' @export
print.rfast99 <- function(x, digits = 4, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y) && ! is.null(x$S)) {
    cat("\nModel runs:", dim(x$y)[1], "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nInteraction indices:\n")
    print(x$I)
    cat("\nTotal order indices:\n")
    print(x$T)
  } else if (! is.null(x$y) && ! is.null(x$mSI)) {
    cat("\nModel runs:", dim(x$y)[1], "\n")
    cat("\n")
    cat("\n==================================")
    cat("\nSensitivity Indices", "\n")
    cat("\nfirst order:", "\n")
    print(round(x$mSI, digits = digits))
    cat("\ninteraction:", "\n")
    print(round(x$iSI, digits = digits))
    cat("\ntotal order:", "\n")
    print(round(x$tSI, digits = digits))
    cat("\n")
    cat("\n=================================")
    cat("\nConvergence Indices", "\n")
    cat("\nfirst order:", "\n")
    print(round(x$mCI, digits = digits))
    cat("\ninteraction:", "\n")
    print(round(x$iCI, digits = digits))
    cat("\ntotal order:", "\n")
    print(round(x$tCI, digits = digits))
  }
  else {
    cat("(empty)\n")
  }
}


#' @rdname rfast99
#' @export
check <- function(x, digits = 4, SI, CI) UseMethod("check")

#' @method check rfast99
#' @export
check.rfast99 <- function(x, digits = 4, SI = 0.01, CI = 0.1){

  if (class(x$mSI)== "array"){
    mSI <- apply(x$mSI, 2, max)
    iSI <- apply(x$iSI, 2, max)
    tSI <- apply(x$tSI, 2, max)
    mCI <- apply(x$mCI, 2, max)
    iCI <- apply(x$iCI, 2, max)
    tCI <- apply(x$tCI, 2, max)
  } else{
    mSI <- x$mSI
    iSI <- x$iSI
    tSI <- x$tSI
    mCI <- x$mCI
    iCI <- x$iCI
    tCI <- x$tCI
  }

  cat("\nSensitivity check ( Index >", SI, ")\n")
  cat("----------------------------------")
  cat("\nfirst order:", names(which(mSI > SI)))
  cat("\ninteraction:", names(which(iSI > SI)))
  cat("\ntotal order:", names(which(tSI > SI)), "\n")
  cat("\n")

  cat("\nConvergence check ( Index >", CI, ")\n")
  cat("----------------------------------")
  cat("\nfirst order:", names(which(mCI > CI)))
  cat("\ninteraction:", names(which(iCI > CI)))
  cat("\ntotal order:", names(which(tCI > CI)))

}

#' @method plot rfast99
#' @importFrom graphics barplot legend lines par abline plot.new
#' @importFrom stats runif fft
#' @export
plot.rfast99 <- function(x, vars = 1, cut.off = F, ...){

  mSI <- x$mSI[,,vars]
  iSI <- x$tSI[,,vars]
  tSI <- x$tSI[,,vars]
  mCI <- x$mCI[,,vars]
  iCI <- x$tCI[,,vars]
  tCI <- x$tCI[,,vars]

  if(is.matrix(mSI)){

    nv <- length(colnames(tSI))+1
    nc <- ceiling(sqrt(nv))
    nr <- ceiling(nv/nc)

    times <- row.names(tSI)

    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc), mar = c(4,2,3,1), oma = c(0,0,2,0))

    for(i in 1:ncol(tSI)){
      plot(times, tSI[,i], ylim = c(0, 1), bty = 'n',
           type = 'l', lwd = 2, xlab = 'time', ylab = '',
           main = colnames(tSI)[i], ...)
      col.transp = adjustcolor('black', alpha.f = 0.4)
      polygon(x = c(times, rev(times)),
              y =c(tSI[,i]-tCI[,i], rev(tSI[,i]+tCI[,i])),
              col = col.transp, border = col.transp)

      col.transp = adjustcolor('red', alpha.f = 0.4)
      lines(times, mSI[,i], ylim = c(0, 1), bty = 'n',
            lwd = 2, col = 'red')
      polygon(x = c(times, rev(times)),
              y =c(mSI[,i]-mCI[,i], rev(mSI[,i]+mCI[,i])),
              col = col.transp, border = col.transp)
      if (is.numeric(cut.off)){
        abline( cut.off, 0, lty = 2)
      }
    }

    variable <- dimnames(x$y)[[4]][vars]
    mtext(variable, NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)

    plot.new()
    legend('top', legend = c('total order', 'first order'), col = c('black','red'),
           lty = 1, lwd = 1, pch = NA, bty = 'n',
           text.col = 'black')
    par(old.par)
  } else {
    D1 <- apply(x$D1, 1, mean)
    V <- apply(x$V, 1, mean)
    Dt <- apply(x$Dt, 1, mean)

    S <- rbind(D1 / V, 1 - Dt / V - D1 / V)
    colnames(S) <- names(x$mSI)
    bar.col <- c("white","grey")
    barplot(S, ylim = c(0,1), col = bar.col)
    legend("topright", c("main effect", "interactions"), fill = bar.col)
  }
}
