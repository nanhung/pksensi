#' Solve PK Model Through MCSim
#'
#' The \code{solve_mcsim} can solve the differential equations of time-dependent quantity/concentration in different tissues/compartments
#' through MCSim.
#'
#' This function allows users to use external data file that assigned in \code{setpoint.name} as parameter matrix.
#' If you want to use it, be sure to define \code{n} and \code{setpoint.name}.
#'
#' @param mName a string giving the name of the model or C file (without extension).
#' @param x a list of storing information in the defined sensitivity function.
#' @param n a numeric to define the sample number.
#' @param dist a vector of distribution names corresponding to \code{<distribution-name>} in MCSim.
#' @param q.arg a list of shape parameters in the sampling distribution (\code{dist}).
#' @param infile.name a character to assign the name of input file.
#' @param setpoint.name a character to assign the name of file for parameter matrix.
#' @param outfile.name a character to assign the name of output file.
#' @param params a character to assign the testing parameters.
#' @param vars a character or a vector to assign the selected output(s).
#' @param time a numeric to define the given time point(s).
#' @param condition a character to set the specific parameter value in the input file.
#' @param rtol an argument passed to the integrator (default 1e-6).
#' @param atol an argument passed to the integrator (default 1e-9).
#' @param generate.infile a logical value to automatically generate the input file, .
#'
#' @importFrom utils write.table
#' @importFrom data.table fread
#'
#' @return The output result is the 4-dimension array with
#' c(model evaluations, replications, time-points, output variables).
#'
#' @examples
#'
#' \dontrun{
#' pbtk1cpt_model()
#' mName <- "pbtk1cpt"
#' compile_model(mName)
#'
#' q <- "qunif"
#' q.arg <- list(list(min = 0.4, max = 1.1),
#'               list(min = 0.1, max = 0.4),
#'               list(min = 1.0, max = 3.0))
#'
#' params <- c("vdist", "ke", "kgutabs")
#'
#' set.seed(1234)
#' x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 20)
#'
#' infile.name <- "example.in"
#' outfile.name <- "example.csv"
#' vars <- "Ccompartment"
#'
#' t <- seq(from = 0.25, to = 12.25, by = 0.5)
#'
#' y <- solve_mcsim(x, mName = mName, infile.name = infile.name,
#'                  setpoint.name = "setpoint.dat",
#'                  outfile.name = outfile.name, params = params, vars = vars, time = t,
#'                  condition = "Agutlument = 10")
#' pksim(y)
#' }
#'
#' @export
#' @describeIn solve_mcsim Numerical analysis for the PK model by MCSim.
solve_mcsim <- function(x, mName,
                        infile.name = NULL,
                        outfile.name = NULL,
                        n = NULL,
                        setpoint.name = NULL,
                        params = NULL,
                        vars  = NULL,
                        time  = NULL,
                        condition  = NULL,
                        generate.infile = T){

  if(is.null(infile.name)) infile.name <- "input.in"
  if(is.null(outfile.name)) outfile.name <- "output.csv"

  if(generate.infile == T){
    generate_infile(infile.name = infile.name,
                    outfile.name = outfile.name,
                    params = params,
                    vars = vars,
                    time = time,
                    condition = condition)
  }

  if(is.null(condition) && is.null(setpoint.name) && is.null(n)){
    stop("Please assign the setpoint.name (parameter matrix defined in input file)")
  }

  if(!is.null(condition)){
    setpoint.data <- "setpoint.dat"
  } else setpoint.data <- setpoint.name

  mcsim. <- paste0("mcsim.", mName)
  #if(file.exists(mcsim.) == F){
  #  stop(paste0("The ", "mcsim.", mName, " doesn't exist."))
  #}

  #
  if (is.numeric(n)){
    n.sample <- n
  } else if (!is.null(x$s)){
    n.sample <- length(x$s)
  }

  if (!is.null(params)){
    n.factors <- length(params)
  } else if (!is.null(x$factors)){
    n.factors <- ifelse(class(x$factors) == "character", length(x$factors), x$factors)
  }

  n.time <- ifelse(is.null(time), 1, length(time))
  n.vars <- length(vars)

  #
  if (is.null(n)){ # Remember to define n if used external parameter matrix
    X <- cbind(1, apply(x$a, 3L, c))
    write.table(X, file=setpoint.data, row.names=F, sep="\t")
  }

  if(file.exists(mcsim.) == T){
    system(paste0("./mcsim.", mName, " ", infile.name))
  } else if (file.exists(paste0(mcsim., ".model.exe")) == T) {
    system(paste0("./mcsim.", mName, ".model.exe ", infile.name))
  }

  if (is.null(n)){rm(X)}

  invisible(gc()); # clean memory

  str <- n.factors + 2
  df <- as.data.frame(data.table::fread(outfile.name, head = T))

  invisible(gc()); # clean memory

  y <- as.matrix(df[,str:ncol(df)]) # output only

  rm(df)
  invisible(gc()); # clean memory

  #  if (!is.null(x)){
  #    n.rep <- x$rep
  #  } else {
  n.rep <- nrow(y) / (n.sample * n.factors)
  #  }

  if (nrow(y) == n.sample){ # For Monte Carlo
    dim <- c(n.sample, 1, n.time, n.vars)
  } else  dim <- c(n.sample * n.factors, n.rep, n.time, n.vars)

  dim(y)<- dim

  invisible(gc()); # clean memory

  dimnames(y)[[3]] <- time
  dimnames(y)[[4]] <- vars

  #file.remove(setpoint.data)

  return(y)
}

#' @export
#' @describeIn solve_mcsim Generate the MCSim input file.
generate_infile <- function(infile.name = NULL,
                            outfile.name = NULL,
                            params, vars, time,
                            condition, rtol = 1e-6, atol = 1e-9,
                            n = NULL, dist = NULL, q.arg = NULL){ # Monte Carlo

  if(is.null(infile.name)) infile.name <- "input.in"
  if(is.null(outfile.name)) outfile.name <- "output.csv"
  setpoint.data <- "setpoint.dat"

  #if(file.exists(paste0(infile.name)) == T){
  #  if(menu(c("Yes", "No"),
  #          title=paste('The "', infile.name, '" is exist. Do you want to replace it?', sep ="")) == 2){
  #    stop()
  #  }
  #}

  cat("#---------------------------------------- \n#",
      " ", infile.name , "\n#",
      " (Created by generate_infile)\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, sep = "")
  cat("Integrate (Lsodes, ", rtol, ", ", atol, " , 1);", "\n\n", file=infile.name, append=TRUE, sep="")
  if(is.null(n)){
    cat("SetPoints (", "\n",
        "\"", outfile.name, "\", \n\"", setpoint.data, "\",\n",
        "0, ", "\n",
        paste(params, collapse = ", "),");\n\n",
        file = infile.name, append = TRUE, sep = "")
  } else {
    cat("MonteCarlo (", "\"", outfile.name, "\"", ",", n , ",", sample(1:99999, 1), ");\n\n",
        file = infile.name, append = TRUE, sep = "")
    for (i in 1 : length(params)){
      cat("Distrib ( ", params[i], ",", dist[i], ",", paste(unlist(q.arg[i]), collapse = ","), ");", "\n",
          file = infile.name, append=TRUE, sep = "")
    }
  }
  cat("\n#---------------------------------------- \n#",
      " Simulation scenario\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, append = TRUE, sep = "")
  cat("Simulation {", "\n\n", file = infile.name, append = TRUE)
  # cat(paste(conditions, collapse=";"), ";", "\n\n", file = infile.name, append=TRUE, sep = "")
  for (i in 1 : length(condition)){
    cat(paste(condition[i], collapse = ";"), ";", "\n", file = infile.name, append = TRUE, sep = "")
  }
  cat("\n", file = infile.name, append = TRUE)
  for (i in 1 : length(vars)) {
    cat("Print (", paste(vars[i], collapse = ", "), ", ", paste(time, collapse=", "), ");\n",
        file = infile.name, append=TRUE, sep = "")
  }
  cat("}", "END.", file = infile.name, append = TRUE)
  message(paste('* Created input file "', infile.name, '".', sep =""))
}
