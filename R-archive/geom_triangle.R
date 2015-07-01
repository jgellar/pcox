geom_triangle <- function(val, limits=NULL, colours=NULL) {
  require(RColorBrewer)
  
  #if (is.null(limits)) limits <- range(val)
  if (is.null(colours)) colours <- rev(brewer.pal(11, "Spectral"))
  #val[val<limits[1]] <- limits[1]
  #val[val>limits[2]] <- limits[2]
  
  geom_tile(aes(fill=val, colour=val)) #+
    #scale_fill_gradientn(name="", limits=limits, colours=colours) +
    #scale_colour_gradientn(name="", limits=limits, colours=colours) +
    #scale_y_continuous(expand = c(0,0)) +
    #scale_x_continuous(expand = c(0,0)) +
    #theme_bw()
}