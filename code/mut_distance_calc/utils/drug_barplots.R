suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# function to create drug barplots
drug_barplots <- function(dat, xlab, ylab, top = 20, fill_var = NULL, title){
  dat <- dat %>%
    dplyr::select(xlab, ylab, fill_var) %>%
    unique() %>%
    filter(get(ylab) != 0) %>%
    arrange(get(ylab)) %>%
    slice_head(n = top) %>%
    as.data.frame()
  
  dat[,xlab] <- factor(dat[,xlab], levels = unique(dat[,xlab]))
  if(!is.null(fill_var)){
    p <- ggplot(dat, aes(x = reorder(get(xlab), -get(ylab)), 
                         y = (-1)*log10(get(ylab)), 
                         fill = get(fill_var))) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("-log10 Adj. P-Value") + 
      scale_fill_manual(name = "Direction", values = c("down" = "forest green", "up" = "red")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
      ggtitle(title)
  } else {
    p <- ggplot(dat, aes(x = get(xlab), 
                         y = (-1)*log10(get(ylab)),
                         fill = (-1)*log10(get(ylab)))) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("-log10 Adj. P-Value") + 
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
      ggtitle(title) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
      guides(fill = "none")
  }
  
  return(p)
}