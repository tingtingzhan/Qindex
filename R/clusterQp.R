

# best explanation!!!
# ~ 1 | c1/c2 
# https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified



#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' Sample \link[stats]{quantile}s in each cluster of observations.
#' 
#' @param formula \link[stats]{formula}, see details in \link[stats]{aggregate.formula}.
#' End user may use following formulas 
#' to specify the response \eqn{y}, cluster(s) \eqn{c}'s and covariate(s) \eqn{x}'s 
#' \describe{
#' \item{`y ~ 1 | c1`}{cluster \eqn{c_1}, without cluster-specific covariate}
#' \item{`y ~ 1 | c1/c2`}{cluster \eqn{c_1}, and cluster \eqn{c_2} nested in \eqn{c_1}, without cluster-specific covariate}
#' \item{`y ~ x1 + x2 | c1`}{cluster \eqn{c_1}, and cluster-specific covariates \eqn{x_1} and \eqn{x_2}}
#' \item{`y ~ x1 + x2 | c1/c2`}{cluster \eqn{c_1}, cluster \eqn{c_2} nested in \eqn{c_1}, and cluster (\eqn{c_2} at least) specific covariates \eqn{x_1} and \eqn{x_2}}
#' \item{`y ~ .  | c1`}{cluster \eqn{c_1}, and all (supposedly cluster-specific) covariates from `data`}
#' \item{`y ~ .  | c1/c2`}{cluster \eqn{c_1}, cluster \eqn{c_2} nested in \eqn{c_1}, and all (supposedly cluster-specific) covariates from `data`}
#' }
#' 
#' @param data \link[base]{data.frame}, 
#' containing 
#' the response \eqn{y}, cluster(s) \eqn{c}'s and covariate(s) \eqn{x}'s 
#' 
# @param FUN \link[base]{character} scalar, name of
# \link[base]{function} for summary statistics, 
# currently only supports \link[stats]{quantile} (`'quantile'` default) 
# \link[np]{npquantile} (`'npquantile'`), or 
# using \link[np]{npquantile} only if sample size is less or equal to 100 (`'npquantile100'`).
#' 
#' @param exclude (optional) \link[base]{character} \link[base]{vector}, 
#' variable(s) excluded from aggregation,
#' e.g., use `exclude = c('z1', 'z2')`
#' to exclude variables \eqn{z_1} and \eqn{z_2}
#' 
#' @param aggregateQp \link[base]{function} used to aggregate the \link[stats]{quantile}s from 
#' cluster \eqn{c_2} (if present). Choices include \link[base]{mean}, \link[base]{max} and \link[base]{min}.
#' 
#' @param from,to,by \link[base]{double} scalars,
#' the starting, end, and increment values
#' to specify a \link[base]{seq}uence of probabilities 
#' \eqn{p = (p_1,\cdots,p_N)'}
#' for the sample \link[stats]{quantile}s \eqn{q = (q_1,\cdots,q_N)'}
#' 
#' @param type \link[base]{integer} scalar, `type` of \link[stats]{quantile} algorithm
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' Function [clusterQp] calculates \eqn{N} sample \link[stats]{quantile}s 
#' in each \link[stats]{aggregate}d cluster of observations.
#' 
#' @returns 
#' Function [clusterQp] returns an \link[stats]{aggregate}d \link[base]{data.frame}.
#' A \link[base]{double} \link[base]{matrix} of \eqn{N} columns is created to store  
#' the sample \link[stats]{quantile}s \eqn{q} of each \link[stats]{aggregate}d cluster.
#' The column names of this \link[stats]{quantile} \link[base]{matrix} are the probabilities \eqn{p}.
#' 
#' @examples 
#' # see ?`Qindex-package` for examples
# @importFrom np npquantile
#' @keywords internal
#' @importFrom stats aggregate quantile update.formula terms.formula
#' @export
clusterQp <- function(
    formula,
    data,
    # FUN = c('npquantile100', 'quantile', 'npquantile'),
    exclude = character(),
    aggregateQp = mean, # max, min
    from = .01, to = .99, by = .01,
    type = 7,
    ...
) {
  
  # if (anyNA(data)) # okay to have NA in `data`
  
  if (length(exclude)) {
    if (!is.character(exclude) || anyNA(exclude)) stop('illegal `exclude`')
    data[exclude] <- NULL
  } #else do nothing
  
  if (!is.symbol(y <- formula[[2L]])) stop('Variable to be aggregated ', sQuote(deparse1(y)), ' must be `symbol`')
  
  if (formula[[3L]][[1L]] != '|') stop('`formula` must have format of y ~ x1 | c1 or y ~ x1 | c1/c2')
  
  # dealing with cluster(s)
  cls_call <- formula[[3L]][[3L]]
  
  # dealing with `x`
  fomx_new <- Reduce(
    f = function(e1, e2) call('-', e1, e2), 
    init = quote(.),
    x = lapply(c(all.vars(y), all.vars(cls_call)), FUN = as.symbol)
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
        data = apply(tmp[[y]], MARGIN = 2L, FUN = aggregateQp, simplify = TRUE),
        dim = c(1L, dim(tmp[[y]])[2L]))
      return(out)
    })
    
    ret <- do.call(what = rbind.data.frame, args = c(tmp2, list(make.row.names = FALSE)))
    
  }
  
  colnames(ret[[y]]) <- probs
  return(ret)
  
}




