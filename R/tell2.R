#' The Decoupling Simulations
#'
#' @description
#' Integrate the decoupling simulations and estimation results
#' to compute the sensitivity measures.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param y a numeric array created from \code{solve_fun} function.
#'
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

#' @rdname tell2
#' @export
tell <- function(x, y = NULL) UseMethod("tell")

#' @method tell rfast99
#' @export
tell.rfast99 <- function(x, y = NULL) {


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
