tidy_index <- function (x, index = "CI") {
  if(!("dplyr" %in% (.packages()))){
    if(!require(dplyr)) install.packages("dplyr") else require(dplyr)
  }
  
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
  names(X) <- c("time", "parameter", "value", "order")
  return(X)
}

ggfast <- function(x, index =  "SI", order = F){
  
  if (index ==  "SI"){
    X <- tidy_index(x, index = index) %>% 
      mutate(category = cut(value, breaks=c(-Inf, 0.01, 0.05, Inf), 
                            labels=c("0 - 0.01","0.01 - 0.05"," > 0.05")))
    
    cols <- c("0 - 0.01" = "grey", "4" = "pink", "0.01 - 0.05" = "pink", " > 0.05" = "red")
  } else if ((index == "CI")) {
    X <- tidy_index(x, index = index) %>% 
      mutate(category = cut(value, breaks=c(-Inf, 0.05, 0.1, Inf), 
                            labels=c("0 - 0.05","0.05 - 0.1"," > 0.1")))
    
    cols <- c("0 - 0.05" = "grey", "4" = "pink", "0.05 - 0.1" = "pink", " > 0.1" = "red")
  }
  
  if (order == F){
    p <- ggplot(X, aes(time, parameter))
  } else if (order == T) {
    p <- ggplot(X, aes(time, reorder(parameter, value)))
  }
  
  p <- p + geom_tile(aes(fill = category), colour = "white") +
    scale_fill_manual(values= cols)+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) + 
    facet_grid(~order) +
    theme(axis.text.x = element_text(size=10, hjust = 1), 
          axis.text.y = element_text(size=10), legend.title=element_blank(),
          legend.position="top")
  
  
  if (index ==  "SI"){
    p + labs(title="Sensitivity index", x="time", y="parameters")
  } else if ((index == "CI")) {
    p + labs(title="Convergence index", x="time", y="parameters")
  }
}
