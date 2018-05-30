#' Create Heatmap to Detect and Check Parameter Sensitivity
#'
#' @description
#' plot the sensitivity (or convergence) index by heatmap with a given result.
#'
#' @param x a list of storing information in the defined sensitivity function.
#' @param fit a vector of interested output index included \code{first order}, \code{interaction}, and \code{total order}.
#' @param vars a logical value or character to specific the display variable in simulation.
#' @param index a character to choose sensitivity index \code{SI} (default) or convergence index \code{CI}.
#' @param order a logical value indicating whether the parameter should reorder by the value.
#' @param level a logical value to use continous or discrete (default) output.
#' @param text a logical value to display the calculated indices in the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom stats reorder time
#' @import dplyr
#' @import ggplot2
#' @export

#' @rdname heat_check
#' @export
heat_check <- function(x, fit = c("first order", "interaction", "total order"),
                       vars = NULL,
                       index = "SI", order = F, level = T, text = F){

  if (index ==  "SI"){
    X <- tidy_index(x, index = index) %>%
      mutate_(level = ~cut(value, breaks=c(-Inf, 0.01, 0.05, Inf),
                           labels=c("0 - 0.01","0.01 - 0.05"," > 0.05")))

    cols <- c("0 - 0.01" = "grey90", "0.01 - 0.05" = "pink1", " > 0.05" = "red")
  } else if ((index == "CI")) {
    X <- tidy_index(x, index = index) %>%
      mutate_(level = ~cut(value, breaks=c(-Inf, 0.05, 0.1, Inf),
                           labels=c("0 - 0.05","0.05 - 0.1"," > 0.1")))

    cols <- c("0 - 0.05" = "grey90", "0.05 - 0.1" = "pink1", " > 0.1" = "red")
  }

  X$variable = factor(X$variable, levels=dimnames(x$y)[[4]])
  X$parameter = factor(X$parameter, levels=rev(x$factors))

  if (is.null(vars)){
    vars <- dimnames(x$y)[[4]]
  } else (vars <- vars)

  X <- X %>% filter(order %in% fit) %>% filter_(~variable %in% vars)

  if (order == F){
    p <- ggplot(X, aes_string("time", "parameter"))
  } else if (order == T) {
    p <- ggplot(X, aes_string("time", "reorder(parameter, value)"))
  }

  if (level == T) {
    p <- p + geom_tile(aes(fill = level), colour = "white") +
      scale_fill_manual(values= cols)
  } else if (level == F) {
    p <- p + geom_tile(aes_string(fill = "value")) +
      scale_fill_gradient(low = "white", high = "red", limits = c(-0.05,1.05))
  }

  p <- p +scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    facet_grid(variable~order) +
    theme(axis.text.x = element_text(size=10, hjust = 1),
          axis.text.y = element_text(size=10), legend.title=element_blank(),
          legend.position="top")

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
