
setOldClass('gam')
setOldClass('gam.prefit')

#' @title Functional Regression Indices
#' 
#' @description
#' Functional regression indices 
#' based on linear or nonlinear functional predictors.
#' 
#' @slot .Data \link[base]{double} \link[base]{vector},
#' functional regression indices, see section **Details**.
#' 
#' @slot formula see explanations in section **Arguments**
#'   
# @slot xgrid strictly increasing \link[base]{double} \link[base]{vector}.
# In package \pkg{Qindex},
# the functional predictor is the \link[stats]{quantile} function, 
# and the input argument for parameter `data` is the returned object of function [clusterQp],
# thus `xgrid`, 
# the common probability grid on which the functional predictor values \eqn{X} are tabulated,
# are the column names of the \link[base]{matrix} \eqn{X}.
#' 
#' @slot gam a \link[mgcv]{gam} object
#' 
#' @slot gpf a `'gam.prefit'` object, which is the returned object 
#' from function \link[mgcv]{gam} with argument `fit = FALSE`
#' 
#' @slot p.value \link[base]{numeric} scalar, 
#' \eqn{p}-value for the test of significance of the functional predictor, 
#' based on slot `@gam`
#' 
#' @slot sign \link[base]{double} scalar of either 1 or -1, 
#' see Step 2 in section **Details**.
#' 
#' @name FRidx
#' @aliases FRidx-class
#' @importFrom methods setClass
#' @export
setClass(Class = 'FRidx', contains = 'numeric', slots = c(
  formula = 'formula',
  gam = 'gam', gpf = 'gam.prefit',
  p.value = 'numeric',
  sign = 'numeric'#, # scalar
  #xgrid = 'numeric' # vector
))


#' @rdname FRidx
#' 
#' @param formula a two-sided \link[stats]{formula} `y ~ X`. 
#' Types of response \eqn{y} may be \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.
#' Functional predictor \eqn{X} is a tabulated \link[base]{double} \link[base]{matrix};
#' the rows of \eqn{X} correspond to the subjects, 
#' while the columns of \eqn{X} correspond to a *common tabulating grid* shared by all subjects.
#' The \link[base]{numeric} values of the grid are in the \link[base]{colnames} of \eqn{X}
#' 
#' @param data \link[base]{data.frame}, e.g., returned from function [clusterQp]
#' 
#' @param sign_prob \link[base]{double} scalar between 0 and 1,
#' probability corresponding to 
#' the selected nearest-even \link[stats]{quantile} in the grid, 
#' which is used to determine the \link[base]{sign} of the (non)linear functional regression indices.
#' Default is `.5`, i.e., the nearest-even \link[stats]{median} of the grid
#' 
#' @param family \link[stats]{family} object, 
#' see function \link[mgcv]{gam}.
#' Default values are
#' \itemize{
#' \item `mgcv::cox.ph()` for \link[survival]{Surv} response \eqn{y};
#' \item `binomial(link = 'logit')` for \link[base]{logical} response \eqn{y};
#' \item `gaussian(link = 'identity')` for \link[base]{double} response \eqn{y}
#' }
#' 
#' @param nonlinear \link[base]{logical} scalar, 
#' whether to use nonlinear or linear functional regression model.
#' Default `FALSE`
#' 
#' @param knot_pct (only when `nonlinear = FALSE`) 
#' positive \link[base]{double} scalar, 
#' percentage of the number of columns of \eqn{X},
#' to be used as `knot_value`.  
#' Default is \eqn{40\%}.
#' If `knot_value` is provided by the end-user, then `knot_pct` is ignored.
#' 
#' @param knot_value (only when `nonlinear = FALSE`) 
#' positive \link[base]{integer} scalar, number of knots 
#' (i.e., parameter `k` in the spline smooth function \link[mgcv]{s})
#' used in \link[mgcv]{gam}.
#' Default is the \link[base]{ceiling} of `knot_pct` of
#' the column dimension of \eqn{X}
#' 
# !!!!! we do not know how to set up number of knots for nonlinear case !!!!!!
#' 
#' @param ... additional parameters, currently not in use.
#' 
#' 
#' @details 
#' 
#' Function [FRidx] calculates 
#' the functional regression indices in the following steps.
#' 
#' \enumerate{
#' 
#' \item Fit a functional regression model (via function \link[mgcv]{gam}) 
#' to the response \eqn{y} 
#' using the functional predictor \eqn{X};
#' 
#' \item Obtain the \link[base]{sign} of the \link[stats]{cor}relation between 
#' \itemize{
#' \item {\eqn{X_{\cdot,j}}, 
#' the user-selected \eqn{j}-th column of functional predictor \eqn{X}.
#' By default, this is the column corresponding to 
#' the \link[stats]{median} of the tabulating grid;}
#' \item `gam(.)$linear.predictors`, from Step 1
#' }
#' 
#' }
#' 
#' *Functional regression indices* (slot `@@.Data`)
#' are the product of 
#' `sign` (from Step 2) and `gam(.)$linear.predictors` (from Step 1).
#' Multiplication by `sign` ensures
#' that the functional regression indices
#' are positively correlated with the user-selected \eqn{X_{\cdot,j}}.
#' 
#' 
#' @returns 
#' 
#' Function [FRidx] returns an \link[base]{S4} class \linkS4class{FRidx} object,
#' the slots of which are described in section **Slots**.
#' 
#' 
#' @references 
#' 
#' Cui, E., Crainiceanu, C. M., & Leroux, A. (2021). 
#' Additive Functional Cox Model. Journal of Computational and Graphical Statistics. 
#' \doi{10.1080/10618600.2020.1853550}
#' 
#' Gellar, J. E., Colantuoni, E., Needham, D. M., & Crainiceanu, C. M. (2015). 
#' Cox regression models with functional covariates for survival data. Statistical Modelling.
#' \doi{10.1177/1471082X14565526}
#' 
#' 
#' @examples 
#' # see ?`Qindex-package`
#' @importFrom stats cor quantile binomial gaussian
#' @importFrom methods new
#' @importFrom mgcv gam cox.ph s ti summary.gam
#' @name FRidx
#' @export
FRidx <- function(
    formula, data,
    sign_prob = .5,
    ...
) {
  
  rhs <- formula[[3L]] # right-hand-side
  X <- data[[rhs]]
  if (!is.symbol(rhs) || !is.matrix(X)) stop('Right-hand-side of `formula` must be a symbol, indicating a matrix column in `data`')
  
  dm <- dim(X)
  xgrid <- as.double(colnames(X))
  if (!is.numeric(xgrid) || anyNA(xgrid) || 
      is.unsorted(xgrid, strictly = TRUE)) {
    stop('`data` needs to be a returned object of function clusterQp()')
  }

  gpf_obj <- FRidx_prefit_(formula = formula, data = data, ...)
  
  gam_obj <- gam(G = gpf_obj, data = data, control = list(keepData = TRUE))
  
  sign_id <- which(xgrid == quantile(xgrid, probs = sign_prob, type = 3L))[1L]
  
  cor_ <- cor(
    x = X[, sign_id], # quantiles of `marker` (at selected quantile) 
    y = gam_obj$linear.predictors, # integration
    use = 'complete.obs'
  ) # scalar
  
  # we use `sign_` to make easy the interpretation
  sign_ <- sign(cor_)
  
  return(new(
    Class = 'FRidx', 
    sign_ * gam_obj$linear.predictors, # functional regression index
    formula = formula,
    #xgrid = xgrid,
    gam = gam_obj, gpf = gpf_obj,
    p.value = summary.gam(gam_obj)$s.table[, 'p-value'],
    sign = sign_
  ))
  
}







#' @rdname FRidx
#' @export
FRidx_prefit_ <- function(
    formula, data,
    family,
    nonlinear = FALSE,
    knot_pct = .4,
    knot_value = ceiling(nxgrid * knot_pct), 
    ...
) {
  
  rhs <- formula[[3L]] # right-hand-side
  X <- data[[rhs]]
  if (!is.symbol(rhs) || !is.matrix(X)) stop('Right-hand-side of `formula` must be a symbol, indicating a matrix column in `data`')
  
  dm <- dim(X)
  xgrid <- as.double(colnames(X))
  nxgrid <- length(xgrid)
  if (!is.numeric(xgrid) || anyNA(xgrid) || 
      is.unsorted(xgrid, strictly = TRUE)) {
    stop('`data` needs to be a returned object of function clusterQp()')
  }
  
  y <- eval(formula[[2L]], envir = data)
  
  xgrid_ <- tcrossprod(rep(1, times = dm[1L]), xgrid)
  
  # for numeric integration of the functional term
  L <- array(1/nxgrid, dim = dm)
  
  trm_ <- if (nonlinear) {
    call(name = 'ti', quote(xgrid_), rhs, by = quote(L), 
         bs = 'cr', # cubic regression spline
         mc = c( # see ?mgcv::ti; which marginals should have centering constraints applied
           FALSE, 
           TRUE
         ), k = 10L) # in Erjia's code
         #), k = 20L) # in our old practice, which is slow
    #call('ti', quote(xgrid_), rhs, by = quote(L), bs = c('cr', 'cr'), mc = c(FALSE, TRUE), k = 20L)
    # we do not know how to specify number of knots for nonlinear case!!!
  } else call(name = 's', quote(xgrid_), by = call('*', quote(L), rhs), bs = 'cr', k = knot_value)
  
  gam_cl <- if (inherits(y, what = 'Surv')) call(
    name = 'gam', 
    formula = call(name = '~', quote(y[,1L]), trm_),
    weights = quote(y[,2L]), 
    family = if (missing(family)) quote(cox.ph()) else substitute(family),
    fit = FALSE, data = quote(data)
  ) else call(
    name = 'gam', 
    formula = call(name = '~', quote(y), trm_),
    family = if (!missing(family)) {
      substitute(family) 
    } else if (is.logical(y) || all(y %in% c(0, 1))) {
      quote(binomial(link = 'logit'))
    } else if (is.numeric(y)) {
      quote(gaussian(link = 'identity'))
    } else stop('not supported yet'), 
    fit = FALSE, data = quote(data)
  )
  
  return(eval(gam_cl)) # class 'gam.prefit'

}














#' @title Show \linkS4class{FRidx} Object
#' 
#' @description
#' Show \linkS4class{FRidx} object.
#' 
#' @param object an \linkS4class{FRidx} object
#' 
#' @returns 
#' The \link[base]{S4} \link[methods]{show} method of \linkS4class{FRidx} object 
#' does not have a returned value.
#' 
#' @keywords internal
#' @importFrom methods show signature
#' @importFrom utils head
#' @export
setMethod(f = show, signature = signature(object = 'FRidx'), definition = function(object) {
  
  cat('p-value from gam: test significance of `marker` as a functional predictor\n')
  print(object@p.value)

  cat('functional regression index\n')
  print(head(c(object)))
    
  return(invisible())
  
})



