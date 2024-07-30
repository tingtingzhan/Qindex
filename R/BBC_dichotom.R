

#' @title Bootstrap-based Optimism Correction for Dichotomization
#' 
#' @description 
#' 
#' Multivariable regression model with bootstrap-based optimism correction on the dichotomized predictors.
#' 
#' @param formula \link[stats]{formula}, e.g., `y~z~x` or `y~1~x`.
#' Types of response \eqn{y} may be \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.
#' Predictors \eqn{x}'s to be dichotomized may be one or more \link[base]{numeric} \link[base]{vector}s and/or one \link[base]{matrix}.
#' Additional predictors \eqn{z}'s, if any, may be of any type.
#' 
#' @param fom \link[stats]{formula}, e.g., `y~z` or `y~1`, for helper functions, with the response \eqn{y} and additional predictors \eqn{z}'s, if any
#' 
#' @param data \link[base]{data.frame}
#'  
#' @param X \link[base]{numeric} \link[base]{matrix} of \eqn{k} columns, 
#' \link[base]{numeric} predictors \eqn{x_1,\cdots,x_k} to be dichotomized
#' 
#' @param X. \link[base]{logical} \link[base]{matrix} \eqn{\tilde{X}} of \eqn{k} columns, 
#' dichotomized predictors \eqn{\tilde{x}_1,\cdots,\tilde{x}_k}
#' 
#' @param R positive \link[base]{integer} scalar, 
#' number of bootstrap replicates \eqn{R}, default `100L`
#'  
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' Function [BBC_dichotom] obtains a multivariable regression model with 
#' bootstrap-based optimism correction on the dichotomized predictors.
#' Specifically,
#' 
#' \enumerate{
#' 
#' \item{Obtain the dichotomizing rules \eqn{\mathbf{\mathcal{D}}} of predictors \eqn{x_1,\cdots,x_k} based on response \eqn{y} (via [m_rpartD]).
#' Multivariable regression (with additional predictors \eqn{z}, if any) 
#' with dichotomized predictors \eqn{\left(\tilde{x}_1,\cdots,\tilde{x}_k\right) = \mathcal{D}\left(x_1,\cdots,x_k\right)} (via helper function [coef_dichotom])
#' is the ***apparent performance***.}
#' 
#' \item{Obtain the bootstrap-based optimism based on \eqn{R} copies of bootstrap samples (via helper function [optimism_dichotom]).
#' The \link[stats]{median} of bootstrap-based optimism over \eqn{R} bootstrap copies
#' is the ***optimism-correction*** of the dichotomized predictors \eqn{\tilde{x}_1,\cdots,\tilde{x}_k}.}
# In future, we may expand the options to include the use of trimmed-mean \link[base]{mean.default}`(, trim)`, etc.
#' 
#' \item{Subtract the optimism-correction (in Step 2) from the apparent performance estimates (in Step 1), 
#' only for \eqn{\tilde{x}_1,\cdots,\tilde{x}_k}. 
#' The apparent performance estimates for additional predictors \eqn{z}'s, if any, are not modified.
#' Neither the variance-covariance (\link[stats]{vcov}) estimates 
#' nor the other regression diagnostics, e.g.,
#' \link[stats]{resid}uals,
#' \link[stats]{logLik}elihood,
#' etc.,
#' of the apparent performance are modified for now.
#' This coefficient-only, partially-modified regression model is  
#' the ***optimism-corrected performance***.
#' }
#' }
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
#' m1 = BBC_dichotom(Surv(futime, death) ~ age + sex + mgus ~ kappa + lambda, 
#'  data = flchain_Circulatory)
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
  apparent_cf <- coef_dichotom(fom = fom, X. = apprent_X, data = data) 
  
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
#' ## Bootstrap-Based Optimism
#' 
#' Helper function [optimism_dichotom] computes the bootstrap-based optimism
#' of the dichotomized predictors. Specifically,
#' 
#' \enumerate{
#' 
#' \item{\eqn{R} copies of bootstrap samples are generated. In the \eqn{j}-th bootstrap sample,
#' \enumerate{
#' \item{obtain the dichotomizing rules \eqn{\mathbf{\mathcal{D}}^{(j)}} of predictors \eqn{x_1^{(j)},\cdots,x_k^{(j)}} based on response \eqn{y^{(j)}} (via [m_rpartD])}
#' \item{multivariable regression (with additional predictors \eqn{z^{(j)}}, if any) coefficient estimates \eqn{\mathbf{\hat{\beta}}^{(j)} = \left(\hat{\beta}_1^{(j)},\cdots,\hat{\beta}_k^{(j)}\right)^t} of 
#' the dichotomized predictors \eqn{\left(\tilde{x}_1^{(j)},\cdots,\tilde{x}_k^{(j)}\right) = \mathcal{D}^{(j)}\left(x_1^{(j)},\cdots,x_k^{(j)}\right)} (via [coef_dichotom]) 
#' are the ***bootstrap performance estimate***.}
#' }
#' }
#' 
#' \item{Dichotomize \eqn{x_1,\cdots,x_k} in the *entire data* using each of the bootstrap rules \eqn{\mathcal{D}^{(1)},\cdots,\mathcal{D}^{(R)}}.
#' Multivariable regression (with additional predictors \eqn{z}, if any) coefficient estimates \eqn{\mathbf{\hat{\beta}}^{[j]} = \left(\hat{\beta}_1^{[j]},\cdots,\hat{\beta}_k^{[j]}\right)^t} of 
#' the dichotomized predictors \eqn{\left(\tilde{x}_1^{[j]},\cdots,\tilde{x}_k^{[j]}\right) = \mathcal{D}^{(j)}\left(x_1,\cdots,x_k\right)} (via [coef_dichotom])
#' are the ***test performance estimate***.}
#' 
#' \item{Difference between the bootstrap and test performance estimates, 
#' an \eqn{R\times k} \link[base]{matrix} of \eqn{\left(\mathbf{\hat{\beta}}^{(1)},\cdots,\mathbf{\hat{\beta}}^{(R)}\right)} minus 
#' another \eqn{R\times k} \link[base]{matrix} of \eqn{\left(\mathbf{\hat{\beta}}^{[1]},\cdots,\mathbf{\hat{\beta}}^{[R]}\right)},
#' are the ***bootstrap-based optimism***.}
#' 
#' }
#' 
#' @section Returns of Helper Functions: 
#' 
#' ## Of helper function [optimism_dichotom]
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
  
  rules <- lapply(bts, FUN = function(i) {
    m_rpartD(y = y[i], X = X[i, , drop = FALSE])
  }) # \eqn{R} copies of dichotomizing rules
  
  boot_dX <- lapply(seq_len(R), FUN = function(i) { # (i = 1L)
    rules[[i]](X[bts[[i]], , drop = FALSE])
  }) # \eqn{R} dichotomizing rules applied to \eqn{R} bootstrap samples, respectively
  
  boot_cf <- lapply(seq_len(R), FUN = function(i) { # (i = 1L)
    coef_dichotom(fom = fom, data = data[bts[[i]], ], X. = boot_dX[[i]])
  }) # regression coef of each `boot_dX` based on each bootstrap sample
  # stopifnot(!anyNA(boot_cf, recursive = TRUE))
  
  test_dX <- lapply(rules, FUN = function(fn) fn(X)) 
  # \eqn{R} dichotomizing rules applied to complete data
  
  test_cf <- lapply(test_dX, FUN = function(i) {
    coef_dichotom(fom = fom, data = data, X. = i)
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
#' ## Multivariable Regression Coefficient Estimates of Dichotomized Predictors \eqn{\tilde{x}}'s
#' 
#' Helper function [coef_dichotom] 
#' fits a multivariable Cox proportional hazards (\link[survival]{coxph}) model for \link[survival]{Surv} response, 
#' logistic (\link[stats]{glm}) regression model for \link[base]{logical} response, 
#' or linear (\link[stats]{lm}) regression model for \link[stats]{gaussian} response,
#' with 
#' the dichotomized predictors \eqn{\tilde{x}_1,\cdots,\tilde{x}_k} as well as
#' the additional predictors \eqn{z}'s.
#' 
#' It is almost inevitable to have duplicates among the dichotomized predictors \eqn{\tilde{x}_1,\cdots,\tilde{x}_k}.
#' In such case, the multivariable model is fitted using the unique \eqn{\tilde{x}}'s.
#' 
#' @section Returns of Helper Functions: 
#' 
#' ## Of helper function [coef_dichotom]
#' 
#' Helper function [coef_dichotom] returns a \link[base]{double} \link[base]{vector} of 
#' the regression \link[stats]{coef}ficients of dichotomized predictors \eqn{\tilde{x}}'s, with \link[base]{attributes}
#' \describe{
#' \item{`attr(,'model')`}{the \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} regression model}
#' }
#' In the case of duplicated \eqn{\tilde{x}}'s, the regression coefficients of the unique \eqn{\tilde{x}}'s are duplicated for those duplicates in \eqn{\tilde{x}}'s.
#' 
#' @importFrom stats lm glm binomial
#' @importFrom survival coxph
#' @importFrom utils tail
#' @rdname BBC_dichotom
#' @export
coef_dichotom <- function(fom, X., data) {
  
  y <- eval(fom[[2L]], envir = data)
  if (anyNA(y)) stop('do not allow missingness in the response, for now')
  
  # identify duplicated columns in `X.`
  dX_orig <- X.
  Xc <- asplit(X., MARGIN = 2L) # list of the columns of `X.`
  dupX <- duplicated.default(Xc)
  if (any(dupX)) {
    ids <- match(x = Xc, table = Xc[!dupX])
    X. <- do.call(cbind, args = Xc[!dupX])
  }
  
  nms <- make.names(dimnames(X.)[[2L]])
  data[nms] <- X.
  
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
