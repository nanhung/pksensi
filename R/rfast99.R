#' Extended Fourier Amplitude Sensitivity Test with Random Phase Shift
#'
#'   rfast99 implements the so-called "extended-FAST" method
#'   (Saltelli et al. 1999). This method allows the estimation of first
#'   order and total Sobol' indices for all the factors.
#'

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

#' @rdname tell.rfast99
#' @export
tell.rfast99 <- function(x, y = NULL, ...) {


  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }

  p <- dim(x$a)[3]
  n <- length(x$s)

  V <- array(numeric(p), dim = c(p, x$rep))
  D1 <- array(numeric(p), dim = c(p, x$rep))
  Dt <- array(numeric(p), dim = c(p, x$rep))

  for (j in 1 : x$rep){
    for (i in 1 : p) {
      l <- seq((i - 1) * n + 1, i * n)
      f <- fft(x$y[l, j], inverse = FALSE)
      Sp <- ( Mod(f[2 : (n / 2)]) / n )^2
      V[i, j] <- 2 * sum(Sp)
      D1[i, j] <- 2 * sum(Sp[(1 : x$M) * x$omega[1]])
      Dt[i, j] <- 2 * sum(Sp[1 : (x$omega[1] / 2)])
    }
  }

  S_original <- apply(D1 / V, 1, mean)
  S_se <- apply(D1 / V, 1, function(x) sqrt(var(x)/length(x)))
  S_min_ci <- apply(D1 / V, 1, quantile, probs= c((1-x$conf)/2))
  S_max_ci <- apply(D1 / V, 1, quantile, probs= c(1-(1-x$conf)/2))
  S <- data.frame(S_original, S_se, S_min_ci, S_max_ci)

  T_original <- apply(1 - Dt / V, 1, mean)
  T_se <- apply(1 - Dt / V, 1, function(x) sqrt(var(x)/length(x)))
  T_min_ci <- apply(1 - Dt / V, 1, quantile, probs= c((1-x$conf)/2))
  T_max_ci <- apply(1 - Dt / V, 1, quantile, probs= c(1-(1-x$conf)/2))
  T <- data.frame(T_original, T_se, T_min_ci, T_max_ci)

  I_original <- apply(1 - Dt / V - D1 / V, 1, mean)
  I_se <- apply(1 - Dt / V - D1 / V, 1, function(x) sqrt(var(x)/length(x)))
  I_min_ci <- apply(1 - Dt / V - D1 / V, 1, quantile, probs= c((1-x$conf)/2))
  I_max_ci <- apply(1 - Dt / V - D1 / V, 1, quantile, probs= c(1-(1-x$conf)/2))
  I <- data.frame(I_original, I_se, I_min_ci, I_max_ci)

  names(I) <- names(S) <- names(T) <- c("original", "std. error", "min. c.i.", "max. c.i.")
  row.names(I) <- row.names(S) <- row.names(T) <- dimnames(x$a)[[3]]

  x$V <- V
  x$D1 <- D1
  x$Dt <- Dt
  x$S <- S
  x$T <- T
  x$I <- I

  assign(id, x, parent.frame())
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

#' @rdname tell2
#' @export
tell2 <- function(x, y){

  id <- deparse(substitute(x))

  if (dim(y)[3] == 1){
    X <- tell.rfast99(x, y[,,1])
    x$mSI <- X$S[,"original"]
    x$iSI <- X$I[,"original"]
    x$tSI <- X$T[,"original"]
    x$mCI <- X$S[,"max. c.i."] - X$S[,"min. c.i."]
    x$iCI <- X$I[,"max. c.i."] - X$I[,"min. c.i."]
    x$tCI <- X$T[,"max. c.i."] - X$T[,"min. c.i."]
  } else {
    for ( i in 1:length(dimnames(y)[[3]])){
      X <- tell.rfast99(x, y[,,i])
      if ( i == 1) { # initialize
        x$mSI <- X$S[,"original"]
        x$iSI <- X$I[,"original"]
        x$tSI <- X$T[,"original"]
        x$mCI <- X$S[,"max. c.i."] - X$S[,"min. c.i."]
        x$iCI <- X$I[,"max. c.i."] - X$I[,"min. c.i."]
        x$tCI <- X$T[,"max. c.i."] - X$T[,"min. c.i."]
      } else { # accumulate
        x$mSI <- rbind(x$mSI, X$S[,"original"])
        x$iSI <- rbind(x$iSI, X$I[,"original"])
        x$tSI <- rbind(x$tSI, X$T[,"original"])
        x$mCI <- rbind(x$mCI, X$S[,"max. c.i."] - X$S[,"min. c.i."])
        x$iCI <- rbind(x$iCI, X$I[,"max. c.i."] - X$I[,"min. c.i."])
        x$tCI <- rbind(x$tCI, X$T[,"max. c.i."] - X$T[,"min. c.i."])
      }
    }
  }

  if (dim(y)[3] == 1){
    names(x$mSI) <- names(x$iSI) <- names(x$tSI) <- names(x$mCI) <- names(x$iCI) <- names(x$tCI) <- dimnames(x$a)[[3]]
  } else {
    rownames(x$mSI) <- rownames(x$iSI) <- rownames(x$tSI) <- rownames(x$mCI) <- rownames(x$iCI) <- rownames(x$tCI) <- dimnames(y)[[3]]
    colnames(x$mSI) <- colnames(x$iSI) <- colnames(x$tSI) <- colnames(x$mCI) <- colnames(x$iCI) <- colnames(x$tCI) <- rownames(x$S)
  }

  x$S<-NULL
  x$I<-NULL
  x$T<-NULL

  assign(id, x, parent.frame())
}

#' @rdname check.rfast99
#' @export
check.rfast99 <- function(x, digits = 4, SI = 0.01, CI = 0.1){

  if (class(x$mSI)== "matrix"){
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

  cat("\nSensitivity check ( index >", SI, ")\n")
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
#' @export
plot.rfast99 <- function(x, cut.off = F, ...){

  nv <- length(colnames(x$tSI))
  nc <- ceiling(sqrt(nv))
  nr <- ceiling(nv/nc)

  times <- row.names(x$tSI)

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(nr, nc), mar = c(4,2,4,1))

  for(i in 1:ncol(x$tSI)){
    plot(times, x$tSI[,i], ylim = c(0, 1), bty = 'n',
         type = 'l', lwd = 2, xlab = 'time', ylab = '',
         main = colnames(x$tSI)[i], ...)
    col.transp = adjustcolor('black', alpha = 0.4)
    polygon(x = c(times, rev(times)),
            y =c(x$tSI[,i]-x$tCI[,i], rev(x$tSI[,i]+x$tCI[,i])),
            col = col.transp, border = col.transp)

    col.transp = adjustcolor('red', alpha = 0.4)
    lines(times, x$mSI[,i], ylim = c(0, 1), bty = 'n',
          lwd = 2, col = 'red')
    polygon(x = c(times, rev(times)),
            y =c(x$mSI[,i]-x$mCI[,i], rev(x$mSI[,i]+x$mCI[,i])),
            col = col.transp, border = col.transp)
    if (is.numeric(cut.off)){
      abline( cut.off, 0, lty = 2)
    }
  }
  legend('top', legend = c('total order', 'first order'), col = c('black','red'),
         lty = 'solid', lwd = 1, pch = NA, bty = 'n',
         text.col = 'black',
         fill = adjustcolor(c('black', 'red'), alpha = 0.4), border = NA, cex = 1.2)
  par(old.par)
}

#' @rdname check
#' @export
check <- function(x, digits = 4, SI, CI) UseMethod("check")

tell <- function(x, y = NULL, ...) UseMethod("tell")
