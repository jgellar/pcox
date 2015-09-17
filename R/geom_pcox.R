#' Plot univariate \code{pcox} coefficients
#' 
#' @param yvar the variable containing the estimate values
#' @param se the variable containing the standard errors. Set to \code{NULL}
#'   to omit the confidence interval
#' @param type type of confidence interval, either a ribbon or line
#' @param line_opts list of additional options for drawing the estimate, to be
#'   passed to \code{geom_line}
#' @param ci_opts list of additional options for drawing the confidence interval,
#'   to be passed to \code{geom_ribbon} or \code{geom_line}
#' @param bw if \code{TRUE}, uses \code{\link[ggplot2]{theme_bw()}}
#' 
#' @examples
#' # Put examples here
#' 

geom_pcox1 <- function(yvar, se=se, type=c("ribbon","line"),
                       line_opts=list(), ci_opts=list(), bw=TRUE) {
  type <- match.arg(type)
  
  # Layer for estimate
  g_est <- list(do.call("geom_line",
                        modifyList(list(mapping=aes_q(y=substitute(yvar)), lwd=2),
                                   line_opts)))
  
  # Layer(s) for confidence interval
  g_ci <- if (is.null(substitute(se))) {
    NULL
  } else if (type=="ribbon") {
    list(do.call("geom_ribbon",
                 modifyList(list(mapping=aes_q(ymin=substitute(yvar-1.96*se),
                                               ymax=substitute(yvar+1.96*se)),
                                 colour="lightgrey", fill="lightgrey"),
                            ci_opts)))
  } else if (type=="line") {
    list(do.call("geom_line",
                 modifyList(list(mapping=aes_q(y=substitute(yvar-1.96*se)),
                                 linetype="dashed"),
                            ci_opts)),
         do.call("geom_line",
                 modifyList(list(mapping=aes_q(y=substitute(yvar+1.96*se)),
                                 linetype="dashed"),
                            ci_opts)))
  } else {
    stop("Unrecognized type for geom_pcox")
  }
  
  # Concatenate and return
  c(g_ci, g_est, if (bw) list(theme_bw()))
}


# geom_myci <- function(yvar, se=se) {
#   list(geom_ribbon(mapping=aes_q(ymin=substitute(yvar-1.96*se),
#                                  ymax=substitute(yvar+1.96*se)),
#                    colour="lightgrey", fill="lightgrey"),
#        geom_line(mapping=aes_q(y=substitute(yvar)), lwd=2),
#        theme_bw())
# }




#' Plot bivariate \code{pcox} coefficients as a heatmap using \code{ggplot2}
#' 
#' This \code{geom} calls \code{\link[ggplot2]{geom_tile}}, and modifies some of
#'   the plot parameters to produce pretty heatmaps of bivariate \code{pcox}
#'   coefficients.
#' 
#' @param z variable to be represented by the color of the heatmap
#' @param colours color palette to use for both the \code{fill} and \code{colour}
#'   aesthetics
#' @param expand.x amount to expand grid in x direction
#' @param expand.y amount to expand grid in y direction
#' @param bw if \code{TRUE}, uses \code{theme_bw()}
#' 
#' @export

geom_pcox2 <- function(z, colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),
                       zlim=NULL, expand.x=c(0,0), expand.y=c(0,0), bw=TRUE,
                       ...) {
  
  
  list(
    ggplot2::geom_tile(ggplot2::aes_q(fill=substitute(z), colour=substitute(z)),
                       ...),
    ggplot2::scale_fill_gradientn(name=NULL, colours=colours, limits=zlim,
                                  oob = scales::squish),
    ggplot2::scale_colour_gradientn(name=NULL, colours=colours, limits=zlim),
    if (bw)
      ggplot2::theme_bw(),
    ggplot2::scale_y_continuous(expand = expand.x),
    ggplot2::scale_x_continuous(expand = expand.y)
  )
}


# fig1 <- ggplot(beta, aes(s,t)) +
#   geom_tile(aes(colour=value, fill=value)) +
#   scale_fill_gradientn(name="", colours=rev(brewer.pal(11,"Spectral"))) +
#   scale_colour_gradientn(name="", colours=rev(brewer.pal(11,"Spectral"))) +
#   theme_bw() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(expand = c(0,0)) +

#   