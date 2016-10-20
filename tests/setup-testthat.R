# Set up variables
set.seed(12354)
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
z <- rnorm(N)
male <- rbinom(N, size = 1, prob=.5)
sind <- seq(0,1,length=J)
X <- genX(N, sind)
K <- 100
Z <- genX(N, seq(0,1,length=K))
L <- 100

# plot function for estimated coefficient surfaces:
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
plotMe <- function(est, lims=range(est$value)) {
  ggplot(est, aes(s, t)) +
    geom_tile(aes(fill=value, colour=value)) +
    theme_bw() +
    scale_fill_gradientn(name="", limits=lims,
      colours=rev(brewer.pal(11,"Spectral"))) +
    scale_colour_gradientn(name="", limits=lims,
      colours=rev(brewer.pal(11,"Spectral"))) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))
}