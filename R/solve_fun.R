solve_fun <- function(x, times = NULL, parameters, initState, dllname,
                      func, jacfunc, initfunc, nout = 1, outnames,
                      model = NULL){
  n <- length(x$s)  
  factors <- ifelse (class(x$factors) == "character", length(x$factors), x$factors) 
  replicate <- x$rep
  out <- ifelse (is.null(times), 1, length(times))
  y <- array(dim = c(n * factors, replicate, out), NA)
  
  if (is.null(model) == TRUE){
    for (k in 1 : dim(y)[3]) { #outputs
      
      # Specific time or variable
      inputs = c(0, times[k])
      
      for (i in 1 : dim(y)[2]) { # replicate
        for (j in 1 : dim(y)[1]) { # Model evaluation
          for (p in 1 : factors) { # input individual factors
            parameters[p] <- x$a[j,i,p]
          }
          
          # Integrate
          tmp <- deSolve::ode(initState, inputs, func = func, parms = parameters, 
                              jacfunc = jacfunc, dllname = dllname, 
                              initfunc = initfunc, nout = nout, outnames = outnames)
          y[j,i,k] <- tmp[2, outnames]
        }
      }
    }
  } else {
    for (i in 1 : dim(y)[2]) { # Replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        
        parameters <- x$a[j,i,]
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
