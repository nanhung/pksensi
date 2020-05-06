#' Check the Parameter Sensitivity
#'
#' Visualize and check the sensitivity (or convergence) measurement with the given result.
#'
#' The convergence of sensitivity indices for each parameter is using the approach proposed
#' by Sarrazin et al. (2016).
#' This method quantitatively assesses the convergence by computing the range of
#' 95\% confidence intervals of the sensitivity indices for each parameter across all data points (time and outputs).
#' Using the global approach based on the heatmap visualization combined with the index "cut-off,"
#' can systematically distinguish between "influential" and "non-influential" parameters (Hsieh et al. 2018).
#'
#' @param x a list of the storing information in the defined sensitivity function.
#' @param order a vector of the interested output index, including \code{first order}, \code{interaction}, and \code{total order}.
#' @param vars a logical value or character to specify the display variable in simulation.
#' @param times a logical value or character to specify the display time in simulation.
#' @param SI.cutoff a value or vector to set the cut-off for sensitivity index. The default is 0.05.
#' @param CI.cutoff a value or vector to set the cut-off for convergence index. The default is 0.05.
#' @param index a character to choose sensitivity index \code{SI} (default) or convergence index \code{CI}.
#' @param level a logical value to use continuous or discrete (default) output.
#' @param text a logical value to display the calculated indices in the plot.
#' @param out a logical value to print the checking result to the console.
#' @param show.all a logical value to show all testing parameters in the heatmap. The default is set to \code{FALSE} to show only the influential parameters.
#' @param ... additional arguments to customize the graphical parameters.
#'
#' @importFrom reshape melt
#' @importFrom magrittr %>%
#' @importFrom stats time
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics barplot legend lines par abline plot.new
#' @importFrom stats runif fft
#' @import ggplot2
#' @import dplyr
#'
#' @return The \code{print} function returns sensitivity and convergence indices
#' with given time-step in the console. The \code{check} method provides the summary of
#' parameter sensitivity and convergence according to the given \code{SI.cutoff} and \code{CI.cutoff}.
#' It can distinguish the influential and non-influential parameter by the providing value
#' of \code{SI.cutoff}. The \code{plot} function can generate the
#' time-course functional outputs of first order and interaction indices for each parameter.
#' The default output is the first model variable. The \code{heat_check} provides a convenient way
#' to visualize and distinguish the influential and non-influential parameter by the setting cut-off.
#' The convergence index can examine the stability of the sensitivity index.
#' To check convergence, be sure to conduct the replication in \code{rfast99}.
#'
#' @references
#' Sarrazin, F., Pianosi, F., & Wagener, T. (2016).
#' Global Sensitivity Analysis of environmental models: Convergence and validation.
#' \emph{Environmental Modelling & Software}, 79, 135–152.
#'
#' Hsieh, N. H., Reisfeld, B., Bois, F. Y., & Chiu, W. A. (2018).
#' Applying a global sensitivity analysis workflow to improve the computational efficiencies
#' in physiologically-based pharmacokinetic modeling.
#' \emph{Frontiers in Pharmacology}, 9, 588.
#'
#' @examples
#' q <- "qunif"
#' q.arg <- list(list(min = 0.6, max = 1),
#'   list(min = 0.5, max = 1.5),
#'   list(min = 0.02, max = 0.3),
#'   list(min = 20, max = 60))
#'
#' params <- c("F","KA","KE","V")
#'
#' set.seed(1234)
#' x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 20)
#'
#' time <- seq(from = 0.25, to = 12.25, by = 0.5)
#' out <- solve_fun(x, model = FFPK, time = time, vars = "output")
#'
#'
#' # Check results of sensitivity measures
#' check(out)
#' plot(out)
#' heat_check(out, show.all = TRUE)
#' heat_check(out, index = "CI")
#'
#' @seealso \code{\link{tell2}}
#'
#' @rdname check
#' @export
check <- function(x, times, vars, SI.cutoff, CI.cutoff, out) UseMethod("check")

#' @method check rfast99
#' @export
check.rfast99 <- function(x, times = NULL, vars = NULL, SI.cutoff = 0.05, CI.cutoff = 0.05, out = TRUE){

  if (is.null(times)) times <- dimnames(x$y)[[3]]
  if (is.null(vars)) vars <- dimnames(x$y)[[4]]

  if (length(times) == 1 && length(vars) == 1) {

    mSI <- x$mSI[times,,vars]
    iSI <- x$iSI[times,,vars]
    tSI <- x$tSI[times,,vars]
    mCI <- x$mCI[times,,vars]
    iCI <- x$iCI[times,,vars]
    tCI <- x$tCI[times,,vars]

  } else if (length(times) == 1 && length(vars) > 1) {

    mSI <- apply(x$mSI[times,,vars], 1, max)
    iSI <- apply(x$iSI[times,,vars], 1, max)
    tSI <- apply(x$tSI[times,,vars], 1, max)
    mCI <- apply(x$mCI[times,,vars], 1, max)
    iCI <- apply(x$iCI[times,,vars], 1, max)
    tCI <- apply(x$tCI[times,,vars], 1, max)

  } else if (length(times) > 1 | length(vars) > 1) {

    mSI <- apply(x$mSI[times,,vars], 2, max)
    iSI <- apply(x$iSI[times,,vars], 2, max)
    tSI <- apply(x$tSI[times,,vars], 2, max)
    mCI <- apply(x$mCI[times,,vars], 2, max)
    iCI <- apply(x$iCI[times,,vars], 2, max)
    tCI <- apply(x$tCI[times,,vars], 2, max)

  } else {

    mSI <- apply(x$mSI, 2, max)
    iSI <- apply(x$iSI, 2, max)
    tSI <- apply(x$tSI, 2, max)
    mCI <- apply(x$mCI, 2, max)
    iCI <- apply(x$iCI, 2, max)
    tCI <- apply(x$tCI, 2, max)

  }

  # else{
  #    mSI <- x$mSI
  #    iSI <- x$iSI
  #    tSI <- x$tSI
  #    mCI <- x$mCI
  #    iCI <- x$iCI
  #    tCI <- x$tCI
  #  }

  if (out) {
    cat("\nSensitivity check ( Index >", SI.cutoff, ")\n")
    cat("----------------------------------")
    cat("\nFirst order:\n", names(which(mSI > SI.cutoff)), "\n")
    cat("\nInteraction:\n", names(which(iSI > SI.cutoff)), "\n")
    cat("\nTotal order:\n", names(which(tSI > SI.cutoff)), "\n")
    cat("\nUnselected factors in total order:\n", names(which(tSI <= SI.cutoff)), "\n")
    cat("\n")

    cat("\nConvergence check ( Index >", CI.cutoff, ")\n")
    cat("----------------------------------")
    cat("\nFirst order:\n", names(which(mCI > CI.cutoff)), "\n")
    cat("\nInteraction:\n", names(which(iCI > CI.cutoff)), "\n")
    cat("\nTotal order:\n", names(which(tCI > CI.cutoff)), "\n")
    cat("\n")
  }

  x <- list(mSI = names(which(mSI > SI.cutoff)), iSI = names(which(iSI > SI.cutoff)), tSI = names(which(tSI > SI.cutoff)),
            mCI = names(which(mCI > CI.cutoff)), iCI = names(which(iCI > CI.cutoff)), tCI = names(which(tCI > CI.cutoff)))
}

#' @rdname check
#' @export
heat_check <- function(x,
                       order = c("first order", "interaction", "total order"),
                       vars = NULL, times = NULL,
                       SI.cutoff = c(0.05, 0.1), CI.cutoff = c(0.05, 0.1),
                       index = "SI", level = T, text = F, show.all = FALSE){

  nSI <- length(SI.cutoff)
  SI.labels<-rep(NA, nSI+1)

  for(i in 1:nSI){
    SI.labels[i+1] <- paste0(SI.cutoff[i]," - ",SI.cutoff[i+1])
  }
  SI.labels[1] <- paste0("0 - ", SI.cutoff[1])
  SI.labels[nSI+1] <- paste0(" > ", SI.cutoff[nSI])

  nCI <- length(CI.cutoff)
  CI.labels<-rep(NA, nCI+1)

  for(i in 1:nCI){
    CI.labels[i+1] <- paste0(CI.cutoff[i]," - ",CI.cutoff[i+1])
  }
  CI.labels[1] <- paste0("0 - ", CI.cutoff[1])
  CI.labels[nCI+1] <- paste0(" > ", CI.cutoff[nCI])

  if (index ==  "SI"){
    X <- tidy_index(x, index = index) %>%
      mutate(level = cut(.data$value, breaks=c(-Inf, paste(SI.cutoff), Inf), labels=SI.labels))

    if (!(show.all == TRUE)) {
      check.out <- check.rfast99(x, vars = vars, SI.cutoff = min(SI.cutoff), out = F)
      X <- X %>% filter(.data$parameter %in% check.out$tSI)
      message(paste0("Display ", length(check.out$tSI), " influential parameters from all ", dim(x$a)[3], " examined parameters."))
    }

  } else if ((index == "CI")) {
    X <- tidy_index(x, index = index) %>%
      mutate(level = cut(.data$value, breaks=c(-Inf, paste(CI.cutoff), Inf), labels=CI.labels))
  }

  colfunc <- colorRampPalette(c("red", "grey90"))

  if (index == "SI"){
    cols <- rev(colfunc(nSI+1))
  } else if (index == "CI") cols <- rev(colfunc(nCI+1))

  X$variable = factor(X$variable, levels=dimnames(x$y)[[4]])
  X$parameter = factor(X$parameter, levels=rev(x$params))

  if (is.null(vars)){
    vars <- dimnames(x$y)[[4]]
  } else (vars <- vars)

  if (is.null(times)){
    times <- dimnames(x$y)[[3]]
  } else (times <- times)


  X <- X %>% filter(order %in% !!(order)) %>% filter(.data$variable %in% vars) %>% filter(time %in% times)

  if(length(times) < 16){
    X$time <- as.factor(X$time)
  }

  #if (order == F){
  p <- ggplot(X, aes_string("time", "parameter"))
  #} else if (order == T) {
  #  p <- ggplot(X, aes_string("time", "reorder(parameter, value)"))
  #}

  if (level == T) {
    p <- p + geom_tile(aes(fill = level), colour = "white") +
      scale_fill_manual(values = cols)
  } else if (level == F) {
    p <- p + geom_tile(aes_string(fill = "value")) +
      scale_fill_gradient(low = "white", high = "red", limits = c(-0.05,1.05))
  }

  if(length(times) < 16){
    p <- p + scale_x_discrete(expand=c(0,0))
  } else p <- p + scale_x_continuous(expand=c(0,0))

  if (length(order) == 1){
    p <- p + scale_y_discrete(expand=c(0,0)) +
      facet_grid(~variable) +
      theme(axis.text.x = element_text(size=10, hjust = 1),
            axis.text.y = element_text(size=10), legend.title=element_blank(),
            legend.position="top")
  } else {
    p <- p + scale_y_discrete(expand=c(0,0)) +
      facet_grid(variable~order) +
      theme(axis.text.x = element_text(size=10, hjust = 1),
            axis.text.y = element_text(size=10), legend.title=element_blank(),
            legend.position="top")
  }

  if (index ==  "SI"){
    p <- p + labs(title="Sensitivity index", x="time", y="parameters")
  } else if ((index == "CI")) {
    p <- p + labs(title="Convergence index", x="time", y="parameters")
  }

  if (text == T){
    if (index ==  "SI"){
      p + geom_text(aes_string(label = "ifelse(value < min(SI.cutoff), '', round(value, 2))"), size = 2.5)
    } else if ((index == "CI")) {
      p + geom_text(aes_string(label = "ifelse(value < min(CI.cutoff), '', round(value, 2))"), size = 2.5)
    }
  } else p
}

tidy_index <- function (x, index = "SI") {

  if(index == "CI") {
    m <- reshape::melt(x$mCI) %>% cbind(order = "first order")
    i <- reshape::melt(x$iCI) %>% cbind(order = "interaction")
    t <- reshape::melt(x$tCI) %>% cbind(order = "total order")
    X <- do.call(rbind, list(m, i, t))
  } else if (index == "SI") {
    m <- reshape::melt(x$mSI) %>% cbind(order = "first order")
    i <- reshape::melt(x$iSI) %>% cbind(order = "interaction")
    t <- reshape::melt(x$tSI) %>% cbind(order = "total order")
    X <- do.call(rbind, list(m, i, t))
  }
  names(X) <- c("time", "parameter", "variable", "value", "order")
  return(X)
}

#' @rdname check
#' @method plot rfast99
#' @export
plot.rfast99 <- function(x, vars = 1, SI.cutoff = 0.1, ...){

  mSI <- x$mSI[,,vars]
  iSI <- x$tSI[,,vars]
  tSI <- x$tSI[,,vars]
  mCI <- x$mCI[,,vars]
  iCI <- x$tCI[,,vars]
  tCI <- x$tCI[,,vars]

  if(is.matrix(mSI)){

    nv <- length(colnames(tSI))+1
    nc <- ceiling(sqrt(nv))
    nr <- ceiling(nv/nc)

    times <- row.names(tSI)

    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc), mar = c(4,2,3,1), oma = c(0,0,2,0))

    for(i in 1:ncol(tSI)){
      plot(times, tSI[,i], ylim = c(0, 1), bty = 'n',
           type = 'l', lwd = 2, xlab = 'time', ylab = '',
           main = colnames(tSI)[i], ...)
      col.transp = adjustcolor('black', alpha.f = 0.4)
      polygon(x = c(times, rev(times)),
              y =c(tSI[,i]-tCI[,i], rev(tSI[,i]+tCI[,i])),
              col = col.transp, border = col.transp)

      col.transp = adjustcolor('red', alpha.f = 0.4)
      lines(times, mSI[,i], ylim = c(0, 1), bty = 'n',
            lwd = 2, col = 'red')
      polygon(x = c(times, rev(times)),
              y =c(mSI[,i]-mCI[,i], rev(mSI[,i]+mCI[,i])),
              col = col.transp, border = col.transp)
      if (is.numeric(SI.cutoff)){
        abline(SI.cutoff, 0, lty = 2)
      }
    }

    if (class(vars) == "character") {vars <- which(dimnames(x$y)[[4]] == vars)}

    variable <- dimnames(x$y)[[4]][vars]
    mtext(variable, NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)

    plot.new()
    legend('top', legend = c('total order', 'first order'), col = c('black','red'),
           lty = 1, lwd = 1, pch = NA, bty = 'n',
           text.col = 'black')
    par(old.par)
  }
}

#' @rdname check
#' @method print rfast99
#' @export
print.rfast99 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y) && ! is.null(x$mSI)) {
    digits <- 4

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

    if (x$rep > 1){ # w/ replication
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
  }
  else {
    cat("(empty)\n")
  }
}
