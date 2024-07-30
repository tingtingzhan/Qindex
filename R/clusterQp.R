
#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' Sample \link[stats]{quantile}s in each cluster of observations.
#' 
#' @param formula \link[stats]{formula}
#' to specify the response \eqn{y}, cluster(s) \eqn{c}'s, 
#' cluster-specific covariate(s) \eqn{x}'s to be retained, and
#' cluster-specific covariate(s) \eqn{z}'s to be removed 
#' from `data`, e.g.,
#' \describe{
#' \item{`y ~ 1 | c1`}{cluster \eqn{c_1}, without cluster-specific covariate}
#' \item{`y ~ 1 | c1/c2`}{cluster \eqn{c_1}, and cluster \eqn{c_2} nested in \eqn{c_1}, without cluster-specific covariate}
#' \item{`y ~ x1 + x2 | c1`}{cluster \eqn{c_1}, and cluster-specific covariates \eqn{x_1} and \eqn{x_2}}
# \item{`y ~ x1 + x2 | c1/c2`}{cluster \eqn{c_1}, cluster \eqn{c_2} nested in \eqn{c_1}, and cluster (\eqn{c_2} at least) specific covariates \eqn{x_1} and \eqn{x_2}}
#' \item{`y ~ . | c1`}{cluster \eqn{c_1}, and all (supposedly cluster-specific) covariates from `data`}
# \item{`y ~ . | c1/c2`}{cluster \eqn{c_1}, cluster \eqn{c_2} nested in \eqn{c_1}, and all (supposedly cluster-specific) covariates from `data`}
#' \item{`y ~ . - z1 - z2 | c1`}{cluster \eqn{c_1}, and all (supposedly cluster-specific) covariates, except for \eqn{z_1} and \eqn{z_2}, from `data`}
#' }
#' 
#' @param data \link[base]{data.frame}
#' 
# @param FUN \link[base]{character} scalar, name of
# \link[base]{function} for summary statistics, 
# currently only supports \link[stats]{quantile} (`'quantile'` default) 
# \link[np]{npquantile} (`'npquantile'`), or 
# using \link[np]{npquantile} only if sample size is less or equal to 100 (`'npquantile100'`).
#' 
#' @param f_sum_ \link[base]{function} to summarize the sample \link[stats]{quantile}s from 
#' lower-level cluster \eqn{c_2} (if present), 
#' such as \link[base]{mean} (default), \link[stats]{median}, \link[base]{max}, \link[base]{min}, etc.
#' 
#' @param from,to,by \link[base]{double} scalars,
#' the starting, end, and increment values
#' of a probability \link[base]{seq}uence
#' \eqn{\mathbf{p} = (p_1,\cdots,p_N)'}
#' *shared for all clusters*,
#' where the cluster-specific sample \link[stats]{quantile}s 
#' \eqn{\mathbf{q} = (q_1,\cdots,q_N)'} of response \eqn{y} are calculated
#' 
#' @param type \link[base]{integer} scalar, see argument `type` of function \link[stats]{quantile}
#' 
#' @param ... additional parameters, currently not in use
#' 
# @details 
# Function [clusterQp] calculates \eqn{N} sample \link[stats]{quantile}s 
# in each \link[stats]{aggregate}d cluster of observations.
#' 
#' @returns 
#' Function [clusterQp] returns an \link[stats]{aggregate}d \link[base]{data.frame},
#' in which 
#' 
#' \itemize{
#' 
#' \item {cluster(s) \eqn{c}'s and cluster-specific covariate(s) \eqn{x}'s 
#' are retained. 
#' \itemize{
#' \item {If the input `formula` takes form of `y ~ . | c1` or `y ~ . - z1 | c1`,
#' then all covariates (except for \eqn{z_1}) are considered cluster-specific thus are all retained;}
#' \item {**Currently**, only the highest cluster \eqn{c_1} is retained. 
#' Sample quantiles from lower-level clusters are summarized using `f_sum_`.}
#' }
#' }
#' 
#' \item {response \eqn{y} is removed; instead,  a \link[base]{double} \link[base]{matrix} of \eqn{N} columns stores
#' the cluster-specific sample \link[stats]{quantile}s \eqn{\mathbf{q}} of the response \eqn{y}.
#' This \link[base]{matrix}
#' \itemize{
#' \item {is named after the \link[base]{parse}d \link[base]{expression} of response \eqn{y} in `formula`;}
#' \item {\link[base]{colnames} are the probabilities \eqn{\mathbf{p}}, for the ease of subsequent programming.}
#' }
#' }
#' 
#' }
#' 
#' 
#' @examples 
#' # see ?`Qindex-package` for examples
# @importFrom np npquantile
#' @keywords internal
#' @importFrom stats aggregate quantile update.formula terms.formula
#' @export
clusterQp <- function(
    formula, data,
    # FUN = c('npquantile100', 'quantile', 'npquantile'),
    f_sum_ = mean, # max, min # old variable name `aggregateQp`
    from = .01, to = .99, by = .01,
    type = 7,
    ...
) {
  
  # if (anyNA(data)) # okay to have NA in `data`
  
  if (!is.symbol(y <- formula[[2L]])) stop('Variable to be aggregated ', sQuote(deparse1(y)), ' must be `symbol`')
  
  if (formula[[3L]][[1L]] != '|') stop('`formula` must have format of y ~ x1 | c1 or y ~ x1 | c1/c2')
  
  # dealing with cluster(s)
  cls_call <- formula[[3L]][[3L]]
  
  # dealing with `x`
  fomx_new <- Reduce(
    f = function(e1, e2) call('-', e1, e2), 
    init = quote(.),
    x = lapply(c(as.character(y), all.vars(cls_call)), FUN = as.symbol)
  )
  fomx <- update.formula(
    old = terms.formula(call('~', formula[[3L]][[2L]]), data = data), 
    new = call('~', fomx_new))
  
  probs <- seq.int(from = from, to = to, by = by)
  
  if (is.symbol(cls_call)) { # 1-level clustering
    
    aggr_fom <- eval(call(name = '~', y,
                     call('+', fomx[[2L]], cls_call)))
    
    # ?stats:::aggregate.formula
    ret <- aggregate(aggr_fom, data = data, FUN = quantile, probs = probs, type = type, names = FALSE)
    
  } else {
    
    if (cls_call[[1L]] != '/') stop('multi-level clustering must be specified as `| c1/c2`')
    if (!is.symbol(cls1 <- cls_call[[2L]]) || !is.symbol(cls2 <- cls_call[[3L]])) stop('2-level clustering must be all-symbol, e.g. `| c1/c2`')
    
    aggr_fom <- eval(call(name = '~', y, Reduce(
      f = function(e1, e2) call('+', e1, e2), 
      x = as.list(cls_call)[-1L], 
      init = fomx[[2L]]
    )))
    
    ret0 <- aggregate(aggr_fom, data = data, FUN = quantile, probs = probs, type = type, names = FALSE)
    
    tmp1 <- split.data.frame(ret0, f = ret0[[cls1]])
    # stopifnot(length(unique(ret0[[cls1]])) == length(tmp))
    
    tmp2 <- lapply(tmp1, FUN = function(tmp) {
      # (tmp = tmp1[[1L]])
      lapply(all.vars(fomx), FUN = function(x) if (!all(duplicated(tmp[[x]])[-1L])) stop(sQuote(x), ' not unique per ', sQuote(cls1)))
      out0 <- tmp
      out0[[cls2]] <- out0[[y]] <- NULL
      out <- unique.data.frame(out0)
      if (nrow(out) != 1) stop('wont happen')
      out[[y]] <- array(
        data = apply(tmp[[y]], MARGIN = 2L, FUN = f_sum_, simplify = TRUE),
        dim = c(1L, dim(tmp[[y]])[2L]))
      return(out)
    })
    
    ret <- do.call(what = rbind.data.frame, args = c(tmp2, list(make.row.names = FALSE)))
    
  }
  
  dimnames(ret[[y]]) <- list(
    if (is.symbol(cls_call)) as.character(eval(cls_call, envir = ret)) else stop('will be quick fix'),
    probs
  )
  
  return(ret)
  
}



# best explanation!!!
# ~ 1 | c1/c2 
# https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified


