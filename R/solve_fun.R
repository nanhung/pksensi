#' Solve PK Model Through \pkg{deSolve} Package or Analytical Function
#'
#' Solve time-dependent quantities/concentrations of different variables in PK model
#' through the imported \code{ode} function in \pkg{deSolve} package.
#' It can also be used to solve the function with analytical solution.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param time a vector to define the given time sequence.
#' @param initParmsfun a character for the given specific initial parameter function.
#' @param initState a vector that define the initial values of state variables for the ODE system.
#' @param dllname a string giving the name of the shared library (without extension)
#' that contains the compiled function.
#' @param func the name of the function in the dynamically loaded shared library.
#' @param initfunc the name of the initialization function (which initialises values of parameters), as provided in dllname.
#' @param outnames the names of output variables calculated in the compiled function \code{func}.
#' @param method method used by integrator (\pkg{deSolve}).
#' @param rtol argument passed to integrator (\pkg{deSolve}).
#' @param atol argument passed to integrator (\pkg{deSolve}).
#' @param model the defined analytical equation with functional output.
#' @param lnparam a logical value that make the statement of the log-transformed parameter (default FALSE).
#' @param vars a character for the selected output.
#' @param tell a logical value to automatically combine the result y to decoupling simulation x.
#' @param ... additional arguments for \code{deSolve::ode} method.
#'
#' @references
#' Soetaert, K. E., Petzoldt, T., & Setzer, R. W. (2010).
#' Solving differential equations in R: package deSolve.
#' \emph{Journal of Statistical Software}, 33(9), 1–25.
#'
#' @examples
#'   q <- "qunif"
#'   q.arg <- list(list(min = 0.6, max = 1.0),
#'    list(min = 0.5, max = 1.5),
#'    list(min = 0.02, max = 0.3),
#'    list(min = 20, max = 60))
#'
#'   params <- c("F","KA","KE","V")
#'
#'   set.seed(1234)
#'   x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 20)
#'
#'   time <- seq(from = 0.25, to = 12.25, by = 0.5)
#'   y <- solve_fun(x, model = FFPK, time = time, vars = "output")
#'
#'   pksim(y) # Visualize uncertainty of model output
#'
#' @seealso \code{\link{pksim}}
#'
#' @export
solve_fun <- function(x, time = NULL, initParmsfun = "initParms", initState, dllname = NULL,
                      func = "derivs", initfunc = "initmod", outnames,
                      method ="lsode", rtol=1e-8, atol=1e-12,
                      model = NULL, lnparam = F, vars = NULL, tell = T, ...){

  message(paste0("Starting time: ", Sys.time()))

  n <- length(x$s)
  no.params <- ifelse (class(x$params) == "character", length(x$params), x$params)
  replicate <- x$rep
  out <- ifelse (is.null(time), 1, length(time))

  if (is.null(vars)) vars <- outnames

  n.vars <- length(vars)
  y <- array(dim = c(n * no.params, replicate, out, n.vars), NA)
  # c(Model Evaluations, replicates, time points, n.vars)

  if (is.null(model)){

    # Specific time or variable
    inputs = c(0, time) # NEED TIME AT ZERO

    params <- rep(0, length(x$params))
    names(params) <- x$params

    for (i in 1 : dim(y)[2]) { # replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        for (p in x$params) {
          params[p] <- ifelse (lnparam == T,  exp(x$a[j,i,p]), x$a[j,i,p])
        }

        if (!is.null(dllname)){
          parms <- do.call(initParmsfun, list(params))
          # Use the initParms function from _inits.R file, if the file had defined
          tmp <- deSolve::ode(initState, inputs, parms = parms, outnames = outnames, nout = length(outnames),
                              dllname = dllname, func = func, initfunc = initfunc, method = method,
                              rtol=rtol, atol=atol, ...)
        } else {
          parms <- params
          tmp <- deSolve::ode(initState, inputs, parms = parms, outnames = outnames, nout = length(outnames),
                              func = func, method = method, rtol=rtol, atol=atol, ...)
        }

        for (l in 1:n.vars){
          for (k in 1 : dim(y)[3]) { # output time
            y[j,i,k,l] <- tmp[k+1, vars[l]] # skip zero
          }
        }
      }
    }
  } else {
    for (i in 1 : dim(y)[2]) { # Replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        if (lnparam == T) { params <- exp(x$a[j,i,])}
        else if (lnparam == F) { params <- x$a[j,i,]}

        if (is.null(time)) tmp <- model(params) else tmp <- model(params, time)

        for (k in 1 : dim(y)[3]) { # Output time
          y[j,i,k,1] <- tmp[k]
        }
      }
    }
  }

  if(dim(y)[3] > 1){
    dimnames(y)[[3]] <- time
    dimnames(y)[[4]] <- vars
  } else { # single time point
    dimnames(y)[[3]] <- list(time)
    dimnames(y)[[4]] <- list(vars)
  }

  if (tell == T){
    tell2(x, y)
  }

  message(paste0("Ending time: ", Sys.time()))

  if (tell == T){
    return(x)
  } else return(y)
}

