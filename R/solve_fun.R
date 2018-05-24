#' Solved PK model with given parameters space
#'
#' @description
#' The \code{solve_fun} solves for the time-dependent quantity/concentration in different tissues
#' through the imported deSolve function. It can also solve the function with analytical function.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param times a vector to define the given time sequence.
#' @param parameters parameters passed to \code{func}.
#' @param getParms a string that call c function to calculate the scaling parameters (default "getParms").
#' @param initState a vector that define the initial values of state variables for the ODE system.
#' @param dllname a string giving the name of the shared library (without extension)
#' that contains the compiled function.
#' @param func the name of the function in the dynamically loaded shared library.
#' @param initfunc the name of the initialisation function (which initialises values of parameters), as provided in dllname.
#' @param outnames the names of output variables calculated in the compiled function \code{func}.
#' @param method method used by integrator (default "lsode").
#' @param rtol argument passed to integrator (default 1e-8).
#' @param atol argument passed to integrator (default 1e-12).
#' @param model the defined analytical equation with functional output.
#' @param lnparam a logical value that make the statement of the log-transformed parameter (default FALSE).
#' @param output a character for the selected output.
#'
#' @rdname solve_fun
#' @export
solve_fun <- function(x, times = NULL, parameters, getParms = "getParms", initState, dllname,
                      func, initfunc, outnames,
                      method ="lsode", rtol=1e-8, atol=1e-12,
                      model = NULL, lnparam = F, output){
  n <- length(x$s)
  no.factors <- ifelse (class(x$factors) == "character", length(x$factors), x$factors)
  replicate <- x$rep
  out <- ifelse (is.null(times), 1, length(times))
  y <- array(dim = c(n * no.factors, replicate, out), NA)

  if (is.null(model) == TRUE){

    # Specific time or variable
    inputs = c(0, times) # NEED TIME AT ZERO

    for (i in 1 : dim(y)[2]) { # replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        for (p in x$factors) {
          parameters[p] <- ifelse (lnparam == T,  exp(x$a[j,i,p]), x$a[j,i,p])
        }
        #parms <- initParms(newParms=parameters)
        parms <- .C(getParms, # "getParms" must actually named in c file
                    as.double(parameters),
                    parms=double(length(parameters)),
                    as.integer(length(parameters)))$parms
        names(parms) <- names(parameters)

        # Integrate
        tmp <- deSolve::ode(initState, inputs, parms = parms, outnames = outnames, nout = length(outnames),
                            dllname = dllname, func = func, initfunc = initfunc, method = method,
                            rtol=rtol, atol=atol)

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



