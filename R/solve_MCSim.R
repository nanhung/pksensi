#' Solved PK model through MCSim
#'
#' @description
#' The \code{solve_MCSim} solves for the differential equations of time-dependent quantity/concentration in different tissues/compartments
#' through MCSim (must be installed).
#' The output result is the 4-dimension array with c(model evaluations, replications, time-points, output variables).
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param mName a string giving the name of the model or C file (without extension).
#' @param infile.name a character to assign the name of in file.
#' @param outfile.name a character to assign the name of out file.
#' @param parameters a vector to assign the testing parameters.
#' @param output a character or a vector to assign the selected output(s).
#' @param time a vector to define the given time sequence.
#' @param condition a character to set the specific parameter value in the in file.
#' @param rtol argument passed to integrator (default 1e-6).
#' @param atol argument passed to integrator (default 1e-9).
#' @param standalone a logic value to create the standalone executable file (default FALSE).
#'
#' @importFrom utils write.table
#' @importFrom data.table fread
#'
#' @rdname solve_MCSim
#' @export
solve_MCSim <- function(x, mName, infile.name, outfile.name,
                        parameters = NULL,
                        output  = NULL,
                        time  = NULL,
                        condition  = NULL){

  if(file.exists(infile.name) == F){
    generate_infile(infile.name = infile.name,
                    outfile.name = outfile.name,
                    parameters = parameters,
                    output = output,
                    time = time,
                    condition = condition)
  }

  mcsim. <- paste0("mcsim.", mName)
  if(file.exists(mcsim.) == F){
    makemcsim(mName)
  }

  #
  n.sample <- length(x$s)
  n.factors <- ifelse(class(x$factors) == "character", length(x$factors), x$factors)
  n.rep <- x$rep
  n.time <- ifelse(is.null(time), 1, length(time))
  n.vars <- length(output)
  dim <- c(n.sample * n.factors, n.rep, n.time, n.vars)

  #
  setpoint.data <- "setpoint.dat"
  X <- cbind(1, apply(x$a, 3L, c))
  write.table(X, file=setpoint.data, row.names=F, sep="\t")
  system(paste0("./mcsim.", mName, " ", infile.name))
  str <- length(x$factors) + 2

  df <- as.data.frame(data.table::fread("setpoint.csv", head = T))

  y <- df[,str:ncol(df)] %>% as.matrix() # output only
  dim(y)<- dim
  dimnames(y)[[3]] <- time
  dimnames(y)[[4]] <- output

  file.remove(setpoint.data)

  return(y)
}

#' @rdname solve_MCSim
#' @export
makemcsim <- function(mName, standalone = F){
  if(standalone == F){
    system(paste0("makemcsim", " ", mName, ".model"))
  } else {system(paste0("makemcsims", " ", mName, ".model"))}
}

#' @rdname solve_MCSim
#' @export
generate_infile <- function(infile.name, outfile.name, parameters, output, time,
                       condition, rtol = 1e-6, atol = 1e-9){

  setpoint.data <- "setpoint.dat"

  cat("#---------------------------------------- \n#",
      " ", infile.name , "\n#",
      " (Created by generate_infile)\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, sep = "")

  cat("Integrate (Lsodes, ", rtol, ", ", atol, " , 1);", "\n\n", file=infile.name,append=TRUE,sep="")
  cat("SetPoints (", "\n",
      "\"", outfile.name, "\", \n\"",setpoint.data,"\",\n",
      "0, ", "\n",
      paste(parameters, collapse=", "),");\n\n",
      file = infile.name, append=TRUE, sep="")

  cat("#---------------------------------------- \n#",
      " Simulation scenario\n#",
      "----------------------------------------", "\n\n",
      file = infile.name, append=TRUE, sep = "")

  cat("Simulation {", "\n\n", file = infile.name, append=TRUE)

  cat(condition, "\n\n", file = infile.name, append=TRUE, sep = "")

  cat("Print (", paste(output, collapse=", "), ", ", paste(time, collapse=", "), ");\n\n",
      "}", "\n\n",
      file = infile.name, append=TRUE, sep="")

  cat("END.", file = infile.name, append=TRUE)
}
