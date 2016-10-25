# plot function for estimated coefficient surfaces:
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

coefplot_2d <- function(coef, lims=range(coef$value)) {
  ggplot(coef, aes(s, t)) +
    geom_tile(aes(fill=value, colour=value)) +
    theme_bw() +
    scale_fill_gradientn(name="", limits=lims,
      colours=rev(brewer.pal(11,"Spectral"))) +
    scale_colour_gradientn(name="", limits=lims,
      colours=rev(brewer.pal(11,"Spectral"))) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))
}
coefplot_1d <- function(coef) {
  ggplot(coef, aes_string(x = colnames(coef)[1], y = "value")) +
    geom_line() +
    geom_ribbon(aes(ymax = value + 2*se, ymin =  value - 2*se), 
      fill=rgb(0,0,0,.1)) +
    theme_bw()
}