#' Create Heatmap to Detect and Check Parameter Sensitivity
#'
#' @description
#' plot the sensitivity (or convergence) index by heatmap with a given result.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param fit a vector of interested output index included \code{first order}, \code{interaction}, and \code{total order}.
#' @param vars a logical value or character to specific the display variable in simulation.
#' @param times a logical value or character to specific the display time in simulation.
#' @param SI.cutoff a vector with two numeric value to specific the cut-off points of sensitivity index in parameter ranking.
#' @param CI.cutoff a vector with two numeric value to specific the cut-off points of convergence index in parameter ranking.
#' @param index a character to choose sensitivity index \code{SI} (default) or convergence index \code{CI}.
#' @param order a logical value indicating whether the parameter should reorder by the value.
#' @param level a logical value to use continous or discrete (default) output.
#' @param text a logical value to display the calculated indices in the plot.
#' @param digits a integer to rounds the values in its first argument.
#' to the specified number of decimal places (default 4).
#' @param SI a numeric value to set the cut-off point for sensitivity index (default 0.1).
#' @param CI a numeric vlaue to set the cut-off point for convergence index (default 0.01).
#'
#' @importFrom reshape melt
#' @importFrom magrittr %>%
#' @importFrom stats reorder time
#' @importFrom grDevices colorRampPalette
#' @import dplyr
#' @import ggplot2
#' @export
#' @details
#' The \code{print} method can used to print sensitivity and convergence indices
#' with given time-step. The \code{check} method provide the summary of
#' parameter sensitivity and convergence. The \code{plot} function provide the
#' time-course functional output of both indices for each parameter.
#' @param ... additional arguments to customize the graphical parameters.

#' @rdname check
#' @export
heat_check <- function(x, fit = c("first order", "total order"),
                       vars = NULL, times = NULL,
                       SI.cutoff = c(0.05, 0.1), CI.cutoff = c(0.05, 0.1),
                       index = "SI", order = F, level = T, text = F){

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
      mutate_(level = ~cut(value, breaks=c(-Inf, paste(SI.cutoff), Inf), labels=SI.labels))

  } else if ((index == "CI")) {
    X <- tidy_index(x, index = index) %>%
      mutate_(level = ~cut(value, breaks=c(-Inf, paste(CI.cutoff), Inf), labels=CI.labels))
  }

  colfunc <- colorRampPalette(c("red", "grey90"))

  if (index == "SI"){
    cols <- rev(colfunc(nSI+1))
  } else if (index == "CI") cols <- rev(colfunc(nCI+1))

  X$variable = factor(X$variable, levels=dimnames(x$y)[[4]])
  X$parameter = factor(X$parameter, levels=rev(x$factors))

  if (is.null(vars)){
    vars <- dimnames(x$y)[[4]]
  } else (vars <- vars)

  if (is.null(times)){
    times <- dimnames(x$y)[[3]]
  } else (times <- times)

  X <- X %>% filter(order %in% fit) %>% filter_(~variable %in% vars) %>% filter(time %in% times)

  if(length(times) < 10){
    X$time <- as.factor(X$time)
  }

  if (order == F){
    p <- ggplot(X, aes_string("time", "parameter"))
  } else if (order == T) {
    p <- ggplot(X, aes_string("time", "reorder(parameter, value)"))
  }

  if (level == T) {
    p <- p + geom_tile(aes(fill = level), colour = "white") +
      scale_fill_manual(values = cols)
  } else if (level == F) {
    p <- p + geom_tile(aes_string(fill = "value")) +
      scale_fill_gradient(low = "white", high = "red", limits = c(-0.05,1.05))
  }

  if(length(times) < 10){
    p <- p + scale_x_discrete(expand=c(0,0))
  } else p <- p + scale_x_continuous(expand=c(0,0))

  if (length(fit) == 1){
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
    p + geom_text(aes_string(label = "ifelse(value < 0.01, '', round(value, 2))"), size = 2.5)
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
#' @export
check <- function(x, times, vars, SI, CI) UseMethod("check")

#' @method check rfast99
#' @export
check.rfast99 <- function(x, times = NULL, vars = NULL, SI = 0.05, CI = 0.05){

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

  cat("\nSensitivity check ( Index >", SI, ")\n")
  cat("----------------------------------")
  cat("\nFirst order:\n", names(which(mSI > SI)), "\n")
  cat("\nInteraction:\n", names(which(iSI > SI)), "\n")
  cat("\nTotal order:\n", names(which(tSI > SI)), "\n")
  cat("\nUnselected factors in total order:\n", names(which(tSI <= SI)), "\n")
  cat("\n")

  cat("\nConvergence check ( Index >", CI, ")\n")
  cat("----------------------------------")
  cat("\nFirst order:\n", names(which(mCI > CI)), "\n")
  cat("\nInteraction:\n", names(which(iCI > CI)), "\n")
  cat("\nTotal order:\n", names(which(tCI > CI)), "\n")
  cat("\n")

}


#' @method plot rfast99
#' @importFrom graphics barplot legend lines par abline plot.new
#' @importFrom stats runif fft
#' @export
plot.rfast99 <- function(x, vars = 1, cut.off = F, ...){

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
      if (is.numeric(cut.off)){
        abline( cut.off, 0, lty = 2)
      }
    }

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
print.rfast99 <- function(x, digits = 4, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y) && ! is.null(x$mSI)) {
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

    if (x$rep > 1){ # without replication
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
