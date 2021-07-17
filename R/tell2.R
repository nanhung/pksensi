#' The Decoupling Simulations
#'
#' Integrate the decoupling simulations (parameter sequences) and estimation results
#' to compute the sensitivity measures.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param y a numeric array generated from the solutions of \code{solve_fun}
#' or \code{solve_mcsim} function.
#'
#' @source
#' This function is based on \code{tell} function in \pkg{sensitivity} package,
#' which is an S3 generic method to estimate sensitivity measures by combining sensitivity object (\code{rfast99}) and external simulation results.
#'
#' @export
tell2 <- function(x, y){

  id <- deparse(substitute(x))
  x$mSI <- x$iSI <- x$tSI <- x$mCI <- x$iCI <- x$tCI <- array(dim = c(dim(y)[3], length(x$params), dim(y)[4]), NA) #c(tim-points, params, variables)
  vars <- dimnames(y)[[4]]

  dimnames(x$mSI)[[2]] <- dimnames(x$iSI)[[2]] <- dimnames(x$tSI)[[2]] <- dimnames(x$mCI)[[2]] <- dimnames(x$iCI)[[2]] <- dimnames(x$tCI)[[2]] <- x$params
  dimnames(x$mSI)[[3]] <- dimnames(x$iSI)[[3]] <- dimnames(x$tSI)[[3]] <- dimnames(x$mCI)[[3]] <- dimnames(x$iCI)[[3]] <- dimnames(x$tCI)[[3]] <- vars
  dimnames(x$mSI)[[1]] <- dimnames(x$iSI)[[1]] <- dimnames(x$tSI)[[1]] <- dimnames(x$mCI)[[1]] <- dimnames(x$iCI)[[1]] <- dimnames(x$tCI)[[1]] <- dimnames(y)[[3]] # time-points

  if (x$replicate == 1) {
    for (k in vars){ # variables
      for ( i in seq_along((dimnames(y)[[3]])) ){  # time-points
        X <- tell.rfast99(x, y[,,i,k])
        x$mSI[i,,k] <- X$S
        x$iSI[i,,k] <- X$I
        x$tSI[i,,k] <- X$T
      }
    }
  } else {
    for (k in vars){ # variables
      for ( i in seq_along((dimnames(y)[[3]])) ){  # time-points
        X <- tell.rfast99(x, y[,,i,k])
        x$mSI[i,,k] <- X$S[,"original"]
        x$iSI[i,,k] <- X$I[,"original"]
        x$tSI[i,,k] <- X$T[,"original"]
        x$mCI[i,,k] <- X$S[,"max. c.i."] - X$S[,"min. c.i."]
        x$iCI[i,,k] <- X$I[,"max. c.i."] - X$I[,"min. c.i."]
        x$tCI[i,,k] <- X$T[,"max. c.i."] - X$T[,"min. c.i."]
      }
    }
  }

  x$y<-y
  x$S<-NULL
  x$I<-NULL
  x$T<-NULL

  assign(id, x, parent.frame())
}

tell.rfast99 <- function(x, y = NULL) {


  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  }

  p <- dim(x$a)[3]
  n <- length(x$s)

  V <- array(numeric(p), dim = c(p, x$replicate))
  D1 <- array(numeric(p), dim = c(p, x$replicate))
  Dt <- array(numeric(p), dim = c(p, x$replicate))

  if (x$replicate == 1){
    for (i in 1 : p) {
      l <- seq((i - 1) * n + 1, i * n)
      f <- fft(x$y[l], inverse = FALSE)
      Sp <- ( Mod(f[2 : (n / 2)]) / n )^2
      V[i] <- 2 * sum(Sp)
      D1[i] <- 2 * sum(Sp[(1 : x$M) * x$omega[1]])
      Dt[i] <- 2 * sum(Sp[1 : (x$omega[1] / 2)])
      S <- apply(D1 / V, 1, mean)
      T <- apply(1 - Dt / V, 1, mean)
      I <- apply(1 - Dt / V - D1 / V, 1, mean)
      names(S) <- names(T) <- names(I) <- dimnames(x$a)[[3]]
    }
  } else {
    for (j in 1 : x$replicate){
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
  }

  x$V <- V
  x$D1 <- D1
  x$Dt <- Dt
  x$S <- S
  x$T <- T
  x$I <- I

  assign(id, x, parent.frame())
}
