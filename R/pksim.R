#' Pharmacokinetic Simulation from Sampling Parameter
#'
#' Pharmacokinetic plot of the output results based on the given parameter (Uncertainty analysis).
#' If the user define the multiple output in model, the generated result will based on
#' first model variable (default). A pharmacokinetic plot with median and the range
#' of min-max, 10%-90%, and 25%-75%.
#'
#' @param y a numeric array created from \code{solve_fun} or \code{solve_mcsim} function.
#' @param vars a logical value or character to specific the display variable in simulation (default 1).
#' @param log a logical value to transform the y-axis to log scale.
#' @param legend a logical value to display the legend in the created plot.
#' @param ... additional arguments to customize the graphical parameters.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics plot polygon axis
#' @importFrom stats quantile
#'
#' @export
pksim <- function(y, vars = 1, log = F, legend = T, ...){

  if (is.null(dim(y))) {
    y <- y$y
  }
  times <- as.numeric(colnames(y[, 1, , vars]))
  if (dim(y)[3] == 1)
    stop("The time point must greater than 1")
  if (dim(y)[2] == 1) {
    quantY <- apply(y[, , , vars], 2, quantile, c(0.5, 0,
                                                  1, 0.1, 0.9, 0.25, 0.75), na.rm = TRUE)
  }
  else quantY <- apply(y[, , , vars], 3, quantile, c(0.5, 0,
                                                     1, 0.1, 0.9, 0.25, 0.75), na.rm = TRUE)
  ytck <- pretty(c(min(quantY, na.rm = TRUE), max(quantY, na.rm = TRUE)))
  if (log == T) {
    quantY[quantY <= 0] <- NA
    quantY <- log(quantY, 10)
    min <- floor(range(quantY)[1])
    max <- ceiling(range(quantY)[2])
    ytck <- seq(min, max, 1)
    natural_ytck <- 10^ytck
  }
  col.transp = adjustcolor("black", alpha.f = 0.2)
  Time <- times
  Output <- quantY[1, ]
  plot(Time, Output, type = "l",
       #xlab = "",  ylab = "",
       ylim = c(min(ytck), max(ytck)), lty = 1,
       las = 1, lwd = 2, col = 1, yaxt = "n", ...)
  if (log == T) {
    axis(2, at = ytck, labels = natural_ytck, las = 2)
  }
  else axis(2, at = ytck, labels = ytck, las = 2)
  if (any(is.na(quantY[2:3, ])) == FALSE) {
    polygon(x = c(times, rev(times)), y = c(quantY[2, ],
                                            quantY[3, seq(from = length(times), to = 1, by = -1)]),
            col = col.transp, lty = 0)
  }
  if (any(is.na(quantY[4:5, ])) == FALSE) {
    polygon(x = c(times, rev(times)), y = c(quantY[4, ],
                                            quantY[5, seq(from = length(times), to = 1, by = -1)]),
            col = col.transp, lty = 0)
  }
  if (any(is.na(quantY[6:7, ])) == FALSE) {
    polygon(x = c(times, rev(times)), y = c(quantY[6, ],
                                            quantY[7, seq(from = length(times), to = 1, by = -1)]),
            col = col.transp, lty = 0)
  }
  if (legend == T) {
    legend("topright", legend = c("Median", "25%-75%",
                                  "10%-90%", "min-max"), col = c("black",
                                                                 "black", "black", "black"), lty = c(1,
                                                                                                     NA, NA, NA), lwd = c(2, NA, NA, NA), pch = NA, bty = "n",
           text.col = "black", fill = adjustcolor(c(NA,
                                                    "black", "grey30", "grey"),
                                                  alpha.f = 0.5), border = NA, cex = 0.7)
  }

}
