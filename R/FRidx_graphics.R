

#' @title Visualize \linkS4class{FRidx} object using R package \pkg{graphics}
#' 
#' @description
#' Create \link[graphics]{persp}ective and \link[graphics]{contour}
#' plots of FR-index integrand using R package \pkg{graphics}.
#' 
#' End users are encouraged to use function [integrandSurface]
#' with \CRANpkg{plotly} work horse.
#' 
#' @param x \linkS4class{FRidx} object
#' 
#' @param n \link[base]{integer} scalar, fineness of visualization,
#' default `501L`. See parameter `n.grid` of function \link[mgcv]{vis.gam}.
#' 
#' @param xlab,ylab \link[base]{character} scalars
#' 
#' @param zlab \link[base]{character} scalar, for function [persp.FRidx]
#' 
#' @param image_col argument `col` of \link[graphics]{image.default}
#' 
#' @param ... ..
#' 
#' @returns
#' Function [persp.FRidx], 
#' a method dispatch of S3 generic \link[graphics]{persp},
#' does not have a return value.
#' 
#' @keywords internal
#' @name FRidx_graphics
#' @importFrom graphics persp
#' @export persp.FRidx
#' @export
persp.FRidx <- function(
    x, 
    n = 31L, 
    xlab = 'Percentages',
    ylab = 'Quantiles',
    zlab = 'Integrand of FR-index',
    ...
) {
  
  z <- z_FRidx(x, n = n)
  # ?graphics:::persp.default
  persp(x = attr(z, which = 'xy', exact = TRUE),
        z = z, 
        xlab = xlab, ylab = ylab, zlab = zlab,
        ...)
  
  return(invisible()) # ?graphics:::persp.default has an invisible return!
}





#' @returns
#' Function [contour.FRidx],
#' a method dispatch of S3 generic \link[graphics]{contour},
#' does not have a return value
#' 
#' @rdname FRidx_graphics
#' @importFrom graphics contour contour.default image.default
#' @importFrom grDevices topo.colors
#' @export contour.FRidx
#' @export
contour.FRidx <- function(
    x, 
    n = 501L,
    image_col = topo.colors(20L),
    xlab = 'Percentages',
    ylab = 'Quantiles',
    ...
) {
  z <- z_FRidx(x, n = n)
  xy <- attr(z, which = 'xy', exact = TRUE)
  
  image.default(
    x = xy, z = z, 
    col = image_col, xlab = xlab, ylab = ylab, ...
  )
  
  contour.default(x = xy, z = z, add = TRUE, ...)
  
  return(invisible())
}



if (FALSE) {
  # mgcv::vis.gam
  debug(mgcv::vis.gam)
  undebug(mgcv::vis.gam)
  #vis.gam(nlfr@gam, view = c("xgrid", "v"), n.grid = 11L, plot.type = "contour", color = "topo")
  mgcv::vis.gam(nlfr@gam, view = c("xgrid", "marker"), n.grid = 101L, plot.type = "contour", color = "topo")
  mgcv::vis.gam(nlfr@gam, view = c("xgrid", "marker"), n.grid = 101L, plot.type = "persp", color = "topo")
  persp(nlfr)
}



# @title (Nonlinear) Functional Regression Weights
# 
# @description
# ..
# 
# 
# @param object an \linkS4class{FRidx} object
# 
# @param n \link[base]{integer} scalar number of points in a grid,
# same as parameter `n` of function \link[mgcv]{plot.gam}, 
# or parameter `n.grid` of function \link[mgcv]{vis.gam}
# 
# @note
# I don't think function [weights_FRidx] should be a S3 method dispatch of S3 generic \link[stats]{weights}.
# 
# @returns 
# For *linear* functional regression model,
# the returned value is a \link[base]{double} \link[base]{vector}.
# 
# For *nonlinear* functional regression model,
# the returned value is a \link[base]{double} \link[base]{matrix}.
# 
# @details
# *Linear functional regression weights*
# are the tabulated weight function on the grid `xarg`.
# These weights are defined as the product of 
# `object@@sign` and the returned value of function \link[mgcv]{predict.gam}.
# 
# @export

#' @importFrom mgcv predict.gam
z_FRidx <- function(
    x, # returned object from [FRidx]
    n = 501L, 
    ...
) {
  
  rhs <- x@formula[[3L]]
  X <- x@gam$data[[rhs]] 
  xgrid <- as.double(colnames(X)) #x@xgrid
  nxgrid <- length(xgrid)
  
  # inspired by ?mgcv::vis.gam
  xy <- list(
    x = seq.int(from = min(xgrid), to = max(xgrid), length.out = n),
    y = seq.int(from = min(X), to = max(X), length.out = n)
  )
  d_surface <- data.frame(
    expand.grid(xy), # span `x` first, then span `y`
    L = 1/nxgrid
  )
  names(d_surface)[1:2] <- c('xgrid_', as.character(rhs))
  
  z0 <- predict.gam(x@gam, newdata = d_surface, se.fit = FALSE, type = 'link')
  z <- array(z0, dim = c(n, n), dimnames = NULL)
  attr(z, which = 'xy') <- xy
  
  return(x@sign * z) # attributes retained
  
}


