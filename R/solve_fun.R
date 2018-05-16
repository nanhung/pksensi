#' @rdname solve_fun
#' @export
solve_fun <- function(x, times = NULL, parameters, initState, dllname,
                      func, initfunc, outnames,
                      method ="lsode", rtol=1e-8, atol=1e-12,
                      model = NULL, lnparam = F, output){
  n <- length(x$s)
  factors <- ifelse (class(x$factors) == "character", length(x$factors), x$factors)
  replicate <- x$rep
  out <- ifelse (is.null(times), 1, length(times))
  y <- array(dim = c(n * factors, replicate, out), NA)

  if (is.null(model) == TRUE){

    # Specific time or variable
    inputs = c(0, times) # NEED TIME AT ZERO

    for (i in 1 : dim(y)[2]) { # replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        for (p in x$factors) {
          parameters[p] <- ifelse (lnparam == T,  exp(x$a[j,i,p]), x$a[j,i,p])
        }

        # Integrate
        tmp <- deSolve::ode(initState, inputs, func = func, parms = parameters,
                            dllname = dllname, method = method, rtol=rtol, atol=atol,
                            initfunc = initfunc, nout = length(outnames), outnames = outnames)

        for (k in 1 : dim(y)[3]) { #outputs
          y[j,i,k] <- tmp[k+1, output] # skip zero
        }
      }
    }
  } else {
    for (i in 1 : dim(y)[2]) { # Replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        if (lnparam == T) { parameters <- exp(x$a[j,i,])}
        else if (lnparam == F) { parameters <- x$a[j,i,]}

        if (is.null(times)) tmp <- model(parameters) else tmp <- model(parameters, times)

        for (k in 1 : dim(y)[3]) { # Output time
          y[j,i,k] <- tmp[k]
        }
      }
    }
  }
  dimnames(y)[[3]]<-times
  return(y)
}


