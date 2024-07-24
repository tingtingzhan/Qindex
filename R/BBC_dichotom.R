

#' @title Bootstrap-based Optimism Correction for Dichotomization
#' 
#' @description 
#' 
#' Functions explained in this documentation are,
#' 
# \describe{
# \item{[BBC_dichotom]}{to obtain a multivariable regression model with bootstrap-based optimism correction on the dichotomized predictors.}
# \item{[optimism_dichotom]}{a helper function to compute the bootstrap-based optimism of the dichotomized predictors.}
# \item{[coef_dichotom]}{a helper function to obtain the estimated multivariable regression coefficients of the dichotomized predictors.}
# }
#' 
#' @param formula \link[stats]{formula} of format `y~z~x` or `y~1~x`.
#' Types of response \eqn{y} may be \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.
#' Predictors \eqn{x}'s to be dichotomized may be one or more \link[base]{numeric} \link[base]{vector}s and/or one \link[base]{matrix}.
#' Additional predictors \eqn{z}'s, if any, may be of any type.
#' 
#' @param fom \link[stats]{formula} of format `y~z` or `y~1`, for helper functions
#' 
#' @param data \link[base]{data.frame}
#'  
#' @param X (for helper function [optimism_dichotom]) 
#' \link[base]{numeric} \link[base]{matrix} of \eqn{k} columns, 
#' \link[base]{numeric} predictors \eqn{x_1,\cdots,x_k} to be dichotomized
#' 
#' @param dX (for helper function [coef_dichotom]) 
#' \link[base]{logical} \link[base]{matrix} of \eqn{k} columns, 
#' a set of \eqn{k} dichotomized predictors
#' 
#' @param R positive \link[base]{integer} scalar, 
#' number of bootstrap replicates \eqn{R}, default `100L`
#'  
#' @param ... additional parameters, currently not in use
#' 
#' 
#' 
#' @details
#' 
#' Function [BBC_dichotom] obtains a multivariable regression model with 
#' bootstrap-based optimism correction on the dichotomized predictors.
#' Specifically,
#' 
#' \enumerate{
#' 
#' \item Dichotomize the \eqn{k} predictors in the *entire data* (using function [m_rpartD]).
#' Fit a regression model to the entire data with the \eqn{k} dichotomized predictors 
#' as well as the additional predictors, if any (using helper function [coef_dichotom]).
#' The estimated regression model is referred to as the *apparent performance*.
#' 
#' \item Obtain the bootstrap-based optimism based on \eqn{R} copies of bootstrap samples,
#' using [optimism_dichotom].
#' Calculate the \link[stats]{median} of bootstrap-based optimism, 
#' specific to each of the dichotomized predictors.
#' In future, we may expand the options to include the use of trimmed-mean 
#' \link[base]{mean.default}`(, trim)`, etc.
#' For now, let's refer to the median optimism as 
#' the *optimism-correction* of the \eqn{k} dichotomized predictors.
#' 
#' }
#' 
#' Subtract the optimism-correction (in Step 2)
#' from the apparent performance estimates (in Step 1), 
#' *only for the \eqn{k} dichotomized predictors*. 
#' The apparent performance estimates for the additional predictors, if any, 
#' are not modified.
#' The variance-covariance (\link[stats]{vcov}) estimates of the apparent performance 
#' is not modified, for now.
#' None of the other regression model diagnostics, such as 
#' \link[stats]{resid}uals,
#' \link[stats]{logLik}elihood,
#' etc.,
#' are modified neither, for now.
#' The coefficient-only, partially-modified regression model is referred to as  
#' the *optimism-corrected performance*.
#' 
#' 

#' 
#' 

#' @returns 
#' 
#' Function [BBC_dichotom] returns a 
#' \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} regression model,
#' with \link[base]{attributes},
#' \describe{
#' \item{`attr(,'optimism')`}{the returned object from [optimism_dichotom]}
#' \item{`attr(,'apparent_cutoff')`}{a \link[base]{double} \link[base]{vector}, 
#' cutoff thresholds for the \eqn{k} predictors in the apparent model}
#' } 
#' 
#' @examples 
#' library(survival)
#' data(flchain, package = 'survival') # see more details from ?survival::flchain
#' head(flchain2 <- within.data.frame(flchain, expr = {
#'   mgus = as.logical(mgus)
#' }))
#' dim(flchain3 <- subset(flchain2, futime > 0)) # required by ?rpart::rpart
#' dim(flchain_Circulatory <- subset(flchain3, chapter == 'Circulatory'))
#' 
#' m1 = BBC_dichotom(Surv(futime, death) ~ age + sex + mgus ~ kappa + lambda, data = flchain_Circulatory)
#' summary(m1)
#' attr(attr(m1, 'optimism'), 'cutoff')
#' attr(m1, 'apparent_cutoff')
#' 
#' @importFrom matrixStats colMedians
#' @importFrom stats model.frame.default na.pass
#' @name BBC_dichotom
#' @export
BBC_dichotom <- function(formula, data, ...) {

  if ((formula[[1L]] != '~') || (length(formula) != 3L)) stop('`formula` must be formula')
  
  fom <- eval(formula[[2L]])
  if (is.symbol(fom) || (fom[[1L]] != '~') || (length(fom) != 3L)) stop('`formula` must be of format `y ~ x1+x2 ~ z1+z2+z3`')
  
  dfom <- call('~', formula[[3L]]) # numeric predictors to be dichotomized
  X <- as.matrix.data.frame(model.frame.default(formula = dfom, data = data, na.action = na.pass))
  if (!is.numeric(X) || !is.matrix(X)) stop(sQuote(deparse1(dfom)), ' must only contain numeric predictors')
  
  y <- eval(fom[[2L]], envir = data)
  
  # Apparent performance 
  m_rule <- m_rpartD(y = y, X = X)
  apprent_X <- m_rule(X)
  apparent_cf <- coef_dichotom(fom = fom, dX = apprent_X, data = data) 
  
  # Bootstrap-based optimism
  optimism <- optimism_dichotom(fom = fom, X = X, data = data, ...) 
  
  medianOpt <- colMedians(optimism, useNames = TRUE, na.rm = TRUE)
  ## later: trimmed-mean ?
  
  ret <- attr(apparent_cf, which = 'model', exact = TRUE)
  ncf <- length(ret$coefficients)
  q <- length(apparent_cf) # number of predictors to be dichotomized
  ret$coefficients[(ncf-q+1L):ncf] <- ret$coefficients[(ncf-q+1L):ncf] - medianOpt
  ## Subtract the mean optimism estimates from the apparent performance 
  ## estimates to obtain the optimism-corrected performance estimates.
  
  # Tingting: we update the `$coefficients` of `ret`
  # so that the Wald-type z-statistics and p-values can be automatically calculated using ?summary
  # We need to update
  # ret$var
  # we still need cov(ret$coefficients, medianOpt)
  # !!!! for now, just leave the variance/covariance as it was !!!
  # end of Tingting
  
  attr(ret, which = 'optimism') <- optimism
  attr(ret, which = 'apparent_cutoff') <- attr(apprent_X, which = 'cutoff', exact = TRUE)
  
  attr(ret, which = 'median_optimism') <- c('Deprecated attribute \'median_optimism\'. Use attr(,\'optimism\') instead')
  attr(ret, which = 'apparent_branch') <- c('Deprecated attribute \'apparent_branch\'. Use attr(,\'apparent_cutoff\') instead')
  
  class(ret) <- c('BBC_dichotom', class(ret))
  return(ret)
  
}











#' @section Details on Helper Functions:
#' 
#' ## [optimism_dichotom]
#' 
#' Function [optimism_dichotom] computes the bootstrap-based optimism
#' of the dichotomized predictors.
#' First, \eqn{R} bootstrap samples are generated,
#' for which the end-user may specify a \link[base]{Random} seed, if needed.
#' Then, 
#' 
#' \enumerate{
#' 
#' \item From each of the \eqn{R} bootstrap samples, 
#' obtain the dichotomizing branches for the \eqn{k} predictors to be dichotomized,
#' using function [m_rpartD]
#' 
#' \item Dichotomize the \eqn{k} predictors in each *bootstrap sample* 
#' using the respective dichotomizing branches from Step 1.
#' The regression coefficients estimate for the \eqn{k} dichotomized predictors
#' (using helper function [coef_dichotom])
#' is referred to as the *bootstrap performance estimate*.
#' 
#' \item Dichotomize the \eqn{k} predictors in the *entire data*
#' using each of the bootstrap dichotomizing branches from Step 1.
#' The regression coefficients estimate for the \eqn{k} dichotomized predictors
#' (using helper function [coef_dichotom])
#' is referred to as the *test performance estimate*.
#' 
#' }
#' 
#' The difference between the bootstrap and test performance estimates, 
#' based on each of the \eqn{R} bootstrap samples,
#' are referred to as the bootstrap-based *optimism* or optimistic bias.
#' 
#' @section Returns of Helper Functions: 
#' 
#' Helper function [optimism_dichotom] returns an \eqn{R\times k} \link[base]{double} \link[base]{matrix} of 
#' bootstrap-based optimism, 
#' with \link[base]{attributes}
#' \describe{
#' \item{`attr(,'cutoff')`}{an \eqn{R\times k} \link[base]{double} \link[base]{matrix}, 
#' the \eqn{R} copies of bootstrap cutoff thresholds for the \eqn{k} predictors.
#' See attribute `'cutoff'` of function [m_rpartD]}
#' }
#'
#' @references 
#' 
#' ## For helper function [optimism_dichotom]
#' 
#' Ewout W. Steyerberg (2009) Clinical Prediction Models.
#' \doi{10.1007/978-0-387-77244-8}
#' 
#' Frank E. Harrell Jr., Kerry L. Lee, Daniel B. Mark. (1996) 
#' Multivariable prognostic models: issues in developing models, evaluating
#' assumptions and adequacy, and measuring and reducing errors.
#' \doi{10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4} 
#' 
#' @rdname BBC_dichotom
#' @export
optimism_dichotom <- function(fom, X, data, R = 100L, ...) {
  
  y <- eval(fom[[2L]], envir = data)
  
  bts <- bootid(n = length(y), R = R) # \eqn{R} copies of 'integer' vectors
  
  parts <- lapply(bts, FUN = function(i) {
    m_rpartD(y = y[i], X = X[i, , drop = FALSE])
  }) # \eqn{R} copies of dichotomizing rules
  
  boot_dX <- lapply(seq_len(R), FUN = function(i) { # (i = 1L)
    parts[[i]](X[bts[[i]], , drop = FALSE])
  }) # \eqn{R} dichotomizing rules applied to \eqn{R} bootstrap samples, respectively
  
  boot_cf <- lapply(seq_len(R), FUN = function(i) { # (i = 1L)
    coef_dichotom(fom = fom, data = data[bts[[i]], ], dX = boot_dX[[i]])
  }) # regression coef of each `boot_dX` based on each bootstrap sample
  # stopifnot(!anyNA(boot_cf, recursive = TRUE))
  
  test_dX <- lapply(parts, FUN = function(fn) fn(X)) 
  # \eqn{R} dichotomizing rules applied to complete data
  
  test_cf <- lapply(test_dX, FUN = function(i) {
    coef_dichotom(fom = fom, data = data, dX = i)
  }) # regression coef of each `test_dX` based on complete data
  # stopifnot(!anyNA(test_cf, recursive = TRUE))
  
  # optimistically biased matrix of coefficients
  # .mapply(FUN = `-`, dots = list(boot_cf, test_cf), MoreArgs = NULL) # slower
  ret <- do.call(rbind, args = boot_cf) - do.call(rbind, args = test_cf) 
  
  attr(ret, which = 'boot_branch') <- c('Deprecated attribute \'boot_branch\'.  Use attr(,\'cutoff\') instead')
  attr(ret, which = 'cutoff') <- do.call(rbind, args = lapply(boot_dX, FUN = attr, which = 'cutoff', exact = TRUE))
  
  return(ret)
  
}









#' @section Details on Helper Functions:
#' 
#' ## [coef_dichotom]
#' 
#' Helper function [coef_dichotom] obtains the 
#' estimated multivariable regression coefficients of the dichotomized predictors.
#' A Cox proportional hazards (\link[survival]{coxph}) regression for \link[survival]{Surv} response, 
#' a logistic (\link[stats]{glm}) regression for \link[base]{logical} response, 
#' or a linear (\link[stats]{lm}) regression for \link[stats]{gaussian} response
#' is performed with 
#' \itemize{
#' \item the dichotomous \link[base]{logical} predictors, given as the (\link[base]{unique}) columns of `dX`, and
#' \item the additional predictors specified in `formula`
#' }
#' 
#' The returned coefficient estimates repeat the corresponding estimates of the unique columns of `dX`.
#' 
#' @section Returns of Helper Functions: 
#' 
#' Helper function [coef_dichotom] returns a \link[base]{double} \link[base]{vector} of the
#' coefficients of the dichotomized predictors, with \link[base]{attributes}
#' \describe{
#' \item{`attr(,'model')`}{the \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} regression model}
#' }
#' 
#'
#' @importFrom stats lm glm binomial
#' @importFrom survival coxph
#' @importFrom utils tail
#' @rdname BBC_dichotom
#' @export
coef_dichotom <- function(fom, dX, data) {
  
  y <- eval(fom[[2L]], envir = data)
  if (anyNA(y)) stop('do not allow missingness in the response, for now')
  
  # identify duplicated columns in `dX`
  dX_orig <- dX
  Xc <- asplit(dX, MARGIN = 2L) # list of the columns of `dX`
  dupX <- duplicated.default(Xc)
  if (any(dupX)) {
    ids <- match(x = Xc, table = Xc[!dupX])
    dX <- do.call(cbind, args = Xc[!dupX])
  }
  
  nms <- make.names(dimnames(dX)[[2L]])
  data[nms] <- dX
  
  fom_orig <- fom # just in case
  for (i in nms) fom[[3]] <- call(name = '+', fom[[3]], as.symbol(i))
  
  suppressWarnings(mod <- if (inherits(y, what = 'Surv')) {
    coxph(formula = fom, data = data)
  } else if (is.logical(y) || all(y %in% c(0, 1))) {
    glm(formula = fom, data = data, family = binomial(link = 'logit'))
  } else lm(formula = fom, data = data))
  
  #if (anyNA(mod$coefficients)) {
  #  print(mod)
  #  warning('still could happen')
  #}
  coef_ <- tail(mod$coefficients, n = length(nms))
  if (any(dupX)) coef_ <- coef_[ids]
  attr(coef_, which = 'model') <- mod
  return(coef_)
  
}




