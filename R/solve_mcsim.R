#' Solve PK Model Through MCSim
#'
#' @description
#' The \code{solve_mcsim} solves for the differential equations of time-dependent quantity/concentration in different tissues/compartments
#' through MCSim.
#' The output result is the 4-dimension array with c(model evaluations, replications, time-points, output variables).
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param n a numeric to define the sample number.
#' @param dist a vector of distribution names corresponding to \code{<distribution-name>} in MCSim.
#' @param q.arg a list of shape parameters in the sampling distribution (\code{dist}).
#' @param mName a string giving the name of the model or C file (without extension).
#' @param infile.name a character to assign the name of input file.
#' @param setpoint.name a character to assign the name of file for parameter matrix.
#' @param outfile.name a character to assign the name of output file.
#' @param parameters a character to assign the testing parameters.
#' @param output a character or a vector to assign the selected output(s).
#' @param time a numeric to define the given time point(s).
#' @param condition a character to set the specific parameter value in the input file.
#' @param rtol argument passed to integrator (default 1e-6).
#' @param atol argument passed to integrator (default 1e-9).
#' @param standalone a logic value to create the standalone executable file (default FALSE).
#'
#' @importFrom utils write.table
#' @importFrom data.table fread
#'
#' @rdname solve_mcsim
#' @export
solve_mcsim <- function(x, mName, infile.name,
                        outfile.name,
                        n = NULL,
                        setpoint.name = NULL,
                        parameters = NULL,
                        output  = NULL,
                        time  = NULL,
                        condition  = NULL){

  if(!is.null(condition)){ # Generate input file if not define condition
    generate_infile(infile.name = infile.name,
                    outfile.name = outfile.name,
                    parameters = parameters,
                    output = output,
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
  if(file.exists(mcsim.) == F){
    makemcsim(paste0(mName, ".model"))
  }

  #
  if (is.numeric(n)){
    n.sample <- n
  } else if (!is.null(x$s)){
    n.sample <- length(x$s)
  }

  if (!is.null(parameters)){
    n.factors <- length(parameters)
  } else if (!is.null(x$factors)){
    n.factors <- ifelse(class(x$factors) == "character", length(x$factors), x$factors)
  }

  n.time <- ifelse(is.null(time), 1, length(time))
  n.vars <- length(output)

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

  y <- as.matrix(df[,str:ncol(df)]) # output only

  #  if (!is.null(x)){
  #    n.rep <- x$rep
  #  } else {
  n.rep <- nrow(y) / (n.sample * n.factors)
  #  }

  if (nrow(df)==n){ # For Monte Carlo
    dim <- c(n.sample, 1, n.time, n.vars)
  } else  dim <- c(n.sample * n.factors, n.rep, n.time, n.vars)

  dim(y)<- dim

  dimnames(y)[[3]] <- time
  dimnames(y)[[4]] <- output

  #file.remove(setpoint.data)

  return(y)
}

#' @rdname solve_mcsim
#' @export
makemcsim <- function(mName, standalone = F){
  if(standalone == F){
    system(paste0("makemcsim", " ", mName, ".model"))
  } else {system(paste0("makemcsims", " ", mName, ".model"))}
}

#' @rdname solve_mcsim
#' @export
generate_infile <- function(infile.name, outfile.name, parameters, output, time,
                            condition, rtol = 1e-6, atol = 1e-9,
                            n = NULL, dist = NULL, q.arg = NULL){ # Monte Carlo

  setpoint.data <- "setpoint.dat"

  cat("#---------------------------------------- \n#",
      " ", infile.name , "\n#",
      " (Created by generate_infile)\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, sep = "")

  cat("Integrate (Lsodes, ", rtol, ", ", atol, " , 1);", "\n\n", file=infile.name,append=TRUE,sep="")

  if(is.null(n)){
    cat("SetPoints (", "\n",
        "\"", outfile.name, "\", \n\"",setpoint.data,"\",\n",
        "0, ", "\n",
        paste(parameters, collapse=", "),");\n\n",
        file = infile.name, append=TRUE, sep="")
  } else {
    cat("MonteCarlo (", "\"", outfile.name, "\"", ",", n , ",", sample(1:99999, 1), ");\n\n",
        file = infile.name, append=TRUE, sep="")
    for (i in 1 : length(parameters)){
      cat("Distrib ( ", parameters[i], ",", dist[i], ",", paste(unlist(q.arg[i]), collapse = ","), ");", "\n",
          file = infile.name, append=TRUE, sep = "")
    }

  }

  cat("\n#---------------------------------------- \n#",
      " Simulation scenario\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, append=TRUE, sep = "")

  cat("Simulation {", "\n\n", file = infile.name, append=TRUE)

  # cat(paste(conditions, collapse=";"), ";", "\n\n", file = infile.name, append=TRUE, sep = "")

  for (i in 1 : length(condition)){
    cat(paste(condition[i], collapse=";"), ";", "\n", file = infile.name, append=TRUE, sep = "")
  }

  cat("\n", file = infile.name, append=TRUE)

  for (i in 1 : length(output)) {
    cat("Print (", paste(output[i], collapse=", "), ", ", paste(time, collapse=", "), ");\n",
        file = infile.name, append=TRUE, sep="")
  }

  cat("}", "END.", file = infile.name, append=TRUE)
}
