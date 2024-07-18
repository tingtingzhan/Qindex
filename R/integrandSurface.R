
#' @title Integrand Surface(s) for One or More \linkS4class{FRidx} Model(s)
#' 
#' @description
#' ..
#' 
#' @param ... one or more \linkS4class{FRidx} objects based on the training set.
#' 
#' @param newdata test \link[base]{data.frame}, with at least 
#' the response \eqn{y^{\text{new}}} and
#' the \link[base]{double} \link[base]{matrix} of 
#' functional predictor values \eqn{X^{\text{new}}}
#' of the test set, tabulated on the same \eqn{x}-grid as the training set.
#' If missing, the training set will be used.
#' 
#' @param n \link[base]{integer} scalar, fineness of visualization,
#' default `501L`. See parameter `n.grid` of function \link[mgcv]{vis.gam}.
#' 
#' @param subj_vis \link[base]{integer} scalar or \link[base]{vector},
#' row indices of `newdata` to be visualized. 
#' Default `1:2`, i.e., the first two test subjects.
#' Use `subj_vis = NULL` to disable visualization of `newdata`.
#' 
#' @param ylim ..
#' 
#' @details 
#' 
#' Function [integrandSurface] ..
#' 
#' 
#' @returns 
#' Function [integrandSurface] returns a pretty \link[plotly]{plotly} object 
#' to showcase the \link[graphics]{persp}ective plot of the
#' ***integrand surface*** of functional regression indices.
#' The ***integrand curves*** of selected test subjects will also be displayed.
#' 
#' 
#' @importFrom mgcv predict.gam
#' @importFrom plotly plot_ly add_paths add_surface layout
#' @importFrom stats asOneSidedFormula predict
#' @export
integrandSurface <- function(
    ...,
    newdata = data,
    # visualization
    n = 501L,
    subj_vis = seq_len(min(50L, .row_names_info(newdata, type = 2L))), 
    ylim = range(X, newX)
) {
  
  xs <- list(...)
  if (!all(vapply(xs, FUN = inherits, what = 'FRidx', FUN.VALUE = NA))) stop('all input needs to be `FRidx`')
  
  rhs_ <- unique(lapply(xs, FUN = function(i) i@formula[[3L]]))
  if (length(rhs_) > 1L) stop('endpoints not same?')
  rhs <- rhs_[[1L]]
  
  data_ <- unique(lapply(xs, FUN = function(i) i@gam$data))
  if (length(data_) > 1L) stop('data not same')
  data <- data_[[1L]]
  
  X <- data[[rhs]]
  
  xgrid <- as.double(colnames(X)) #fr@xgrid
  nxgrid <- length(xgrid)
  
  newX <- newdata[[rhs]]
  if (!is.matrix(newX)) stop('`newdata` does not contain a matrix column of functional predictor values')
  new_xgrid <- as.double(colnames(newX))
  if (!all.equal.numeric(new_xgrid, xgrid)) stop('grid of training and test data must be exactly the same')
  
  # plot!!
  # *surface* based on training model
  x_ <- seq.int(from = min(xgrid), to = max(xgrid), length.out = n)
  y_ <- seq.int(from = ylim[1L], to = ylim[2L], length.out = n)
  d_surface <- data.frame(
    expand.grid(xgrid_ = x_, y_), # span `x_` first, then span `y_`
    L = 1/nxgrid
  )
  names(d_surface)[2] <- c(as.character(rhs))
  zs <- lapply(xs, FUN = function(i) {
    y0 <- i@sign * predict.gam(i@gam, newdata = d_surface, se.fit = FALSE, type = 'link')
    t.default(array(y0, dim = c(n, n), dimnames = NULL))
    # ?base::t.default important here!!!
    # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
  })
  
  zmin <- min(unlist(zs))
  zmax <- max(unlist(zs))
  
  col_scale <- list(
    c(0, 1), 
    # c('lightyellow', 'lightpink') # nice
    c('beige', 'lightpink') # nice
    # c('white', 'deeppink') # not good!
    # c('white', 'magenta') # not good!
    # c('white', 'lightgreen') # nice
    # c('white', 'lightgoldenrod') # my R do not recognize
    # c('white', 'lightslateblue') # my R do not recognize
    # c('white', 'yellow') # nice
  )
  
  p <- plot_ly(x = x_, y = y_)
  for (z_ in zs) {
    p <- add_surface(
      p = p, 
      z = z_, 
      cmin = zmin, cmax = zmax, 
      colorscale = col_scale, 
      showscale = FALSE)
  }
  p <- layout(p = p, scene = list(
    xaxis = list(title = 'Probability', tickformat = '.0%', color = 'dodgerblue'), 
    yaxis = list(title = 'Quantile', color = 'deeppink'),
    zaxis = list(title = 'Integrand of FR-Index', color = 'darkolivegreen')
  ))
  
  if (!length(subj_vis)) return(p)
  
  if (!is.integer(subj_vis) || anyNA(subj_vis) || any(subj_vis > nrow(newX))) stop('illegal `subj_vis`')
  
  d_subj <- data.frame(
    xgrid_ = xgrid,
    y = c(t.default(newX[subj_vis, , drop = FALSE])),
    id = rep(subj_vis, each = nxgrid),
    L = 1/nxgrid
  )
  names(d_subj)[2] <- as.character(rhs)
  z_subj <- lapply(xs, FUN = function(i) {
    i@sign * predict.gam(i@gam, newdata = d_subj, se.fit = FALSE, type = 'link')
  })
  
  for (i in seq_along(xs)) {
    p <- add_paths(
      p = p, data = d_subj, 
      x = ~ xgrid_, y = asOneSidedFormula(rhs), 
      z = z_subj[[i]], 
      name = ~id, color = ~id,
      showlegend = FALSE, # (length(subj_vis) <= 10L),
      line = list(
        width = if (length(subj_vis) <= 5L) 5 else 2
      ))
  }
  
  #p <- layout(p = p, legend = list(
  #  title = list(text = if (identical(newdata, data)) 'Training Subj' else 'Test Subj'))
  #)
  return(p)
  
}






