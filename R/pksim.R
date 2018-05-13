#' @importFrom grDevices adjustcolor
#' @importFrom graphics plot polygon
#' @importFrom stats quantile
#'
#' @rdname pksim
#' @export
pksim <- function(y, log = F, legend = T, ...){
  times <- as.numeric(colnames(y[,1,]))
  quantY <- apply(y, 3, quantile, c(0.50, 0, 1, 0.1, 0.9, 0.25,0.75), na.rm=TRUE)
  ytck <- pretty(c(min(quantY,na.rm=TRUE),max(quantY,na.rm=TRUE)))

  if (log == T){
    quantY[quantY <= 0] <- NA # prevent numerical error
    quantY <- log(quantY)
    ytck <- pretty(c(min(quantY, na.rm=TRUE),max(quantY, na.rm=TRUE)))
  }

  col.transp = adjustcolor('black', alpha.f = 0.2)
  plot(times, quantY[1,], type="l", xlab="", ylab="",
       ylim=c(min(ytck),max(ytck)), lty=1, las=1, lwd=2, col = 1, ...)
  if (any(is.na(quantY[2:3,])) == FALSE)
  {
    polygon(x = c(times, rev(times)), y = c(quantY[2,],quantY[3,seq(from=length(times),to=1,by=-1)]),col=col.transp,lty=0)
  }
  if (any(is.na(quantY[4:5,])) == FALSE)
  {
    polygon(x = c(times, rev(times)), y = c(quantY[4,],quantY[5,seq(from=length(times),to=1,by=-1)]),col=col.transp,lty=0)
  }
  if (any(is.na(quantY[6:7,])) == FALSE)
  {
    polygon(x = c(times, rev(times)), y = c(quantY[6,],quantY[7,seq(from=length(times),to=1,by=-1)]),col=col.transp,lty=0)
  }

  if (legend == T){
    legend('topright', # inset=c(0,1), xpd=TRUE, horiz=TRUE,
           legend = c('Median', '25%-75%', '10%-90%', 'min-max'),
           col = c('black','black','black','black'),
           lty = c(1,NA,NA,NA),
           lwd = c(2,NA,NA,NA),
           pch = NA, bty = 'n', text.col = 'black',
           fill = adjustcolor(c(NA, 'black', 'grey30','grey'),
                              alpha.f = 0.5), border = NA, cex = 0.7)
  }
}
