


#' @title Dichotomize via Recursive Partitioning
#' 
#' @description
#' Dichotomize one or more predictors of
#' a \link[survival]{Surv}, a \link[base]{logical}, or a \link[base]{double} response,
#' using recursive partitioning and regression tree \link[rpart]{rpart}.
#' 
#' @param y a \link[survival]{Surv} object, 
#' a \link[base]{logical} \link[base]{vector}, 
#' or a \link[base]{double} \link[base]{vector}, the response \eqn{y}
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, one predictor \eqn{x}
#' 
#' @param X \link[base]{numeric} \link[base]{matrix}, 
#' a set of predictors.
#' Each column of \eqn{X} is one predictor.
#' 
#' @param check_degeneracy \link[base]{logical} scalar, whether to allow the 
#' dichotomized value to be all-`FALSE` or all-`TRUE` (i.e., degenerate) 
#' for any one of the predictors.
#' Default `TRUE` to produce a \link[base]{warning} message for degeneracy.
#' 
#' @param cp \link[base]{double} scalar, complexity parameter, see \link[rpart]{rpart.control}.
#' Default `.Machine$double.eps`, so that a split is enforced 
#' no matter how small improvement in overall \eqn{R^2} is
#' 
#' @param maxdepth positive \link[base]{integer} scalar, maximum depth of any node, see \link[rpart]{rpart.control}.
#' Default `2L`, because only the first node is needed
#' 
#' @param ... additional parameters of \link[rpart]{rpart} and/or \link[rpart]{rpart.control}
#' 
#' 
#' 
#' @details
#' 
#' ## Dichotomize Single Predictor
#' 
#' Function [rpartD] dichotomizes one predictor in the following steps, 
#' 
#' \enumerate{
#' 
#' \item {Recursive partitioning and regression tree \link[rpart]{rpart} analysis is 
#' performed for the response \eqn{y} and the predictor \eqn{x}.}
#' 
#' \item {The \link[rpart]{labels.rpart} of the first node of 
#' the \link[rpart]{rpart} tree
#' is considered as the dichotomizing rule of the \link[base]{double} predictor \eqn{x}.
#' The term *dichotomizing rule* indicates the combination of an inequality sign
#' (\link[base]{>}, \link[base]{>=}, \link[base]{<} and \link[base]{<=}) 
#' and a \link[base]{double} cutoff threshold \eqn{a}}
#' 
#' \item {The dichotomizing rule from Step 2 is further processed, such that
#' \itemize{
#' \item {\eqn{<a} is regarded as \eqn{\geq a}}
#' \item {\eqn{\leq a} is regarded as \eqn{>a}}
#' \item {\eqn{> a} and \eqn{\geq a} are regarded as is.}
#' }
#' This step is necessary for a narrative of 
#' *greater than* or *greater than or equal to* 
#' the threshold \eqn{a}.}
#' 
#' \item {A \link[base]{warning} message is produced, 
#' if the dichotomizing rule, applied to a new \link[base]{double} predictor `newx`, creates 
#' an all-`TRUE` or all-`FALSE` result.
#' We do not make the algorithm \link[base]{stop}, 
#' as most regression models in R are capable of handling 
#' an all-`TRUE` or all-`FALSE` predictor,
#' by returning a `NA_real_` regression coefficient estimate.
#' }
#' 
#' }
#' 
#' 
#' @returns 
#' 
#' ## Dichotomize Single Predictor
#' 
#' Function [rpartD] returns a \link[base]{function}, 
#' with a \link[base]{double} \link[base]{vector} parameter `newx`.
#' The returned value of `rpartD(y,x)(newx)` is a 
#' \link[base]{logical} \link[base]{vector}
#' with \link[base]{attributes}
#' \describe{
#' \item{`attr(,'cutoff')`}{\link[base]{double} scalar, the cutoff value for `newx`}
#' }
#' 
#' 

#' @note
#' In future \link[base]{integer} and \link[base]{factor} predictors will be supported.
#' 
#' @examples
#' ## Dichotomize Single Predictor
#' data(cu.summary, package = 'rpart') # see more details from ?rpart::cu.summary
#' with(cu.summary, rpartD(y = Price, x = Mileage, check_degeneracy = FALSE))
#' (foo = with(cu.summary, rpartD(y = Price, x = Mileage)))
#' foo(rnorm(10, mean = 24.5))
#' @keywords internal
#' @importFrom rpart rpart
#' @name rpartD
#' @export
rpartD <- function(
    y, x, check_degeneracy = TRUE,
    cp = .Machine$double.eps, # to force a split even if the overall lack of fit is not decreased
    maxdepth = 2L, # only the first node is needed
    ...
) {
  tree <- rpart(formula = y ~ x, data = data.frame(y = y, x = x), cp = cp, maxdepth = maxdepth, ...)
  if (!length(tree$splits)) stop('we must force a split')
  
  labs <- labels(tree) # ?rpart:::labels.rpart
  node1 <- str2lang(labs[2L]) # first node!!!
  
  if (node1[[1L]] == '<=') {
    node1[[1L]] <- quote(`>`)
  } else if (node1[[1L]] == '<') {
    node1[[1L]] <- quote(`>=`)
  } # else if (node1[[1L]] is '>' or '>=')  do nothing
  
  #if (!identical(node1[[2L]], quote(x))) stop('rpart package updated?')
  node1[[2L]] <- quote(newx)
  
  node1[[3L]] <- tree$splits[1L, 4L] # threshold, in case `labels` are truncated due to `digits`
  
  fun <- alist(newx = )
  fun[[2L]] <- if (check_degeneracy) call(
    name = '{',
    call('<-', quote(ret), call('(', node1)),
    quote(if (all(ret, na.rm = TRUE) || !any(ret, na.rm = TRUE)) warning('Dichotomized value is all-0 or all-1')),
    call('<-', call('attr', quote(ret), which = 'cutoff'), node1[[3L]]),
    call('<-', call('attr', quote(ret), which = 'text'), deparse1(node1[c(1L, 3L)])),
    quote(return(ret))
  ) else node1
  return(as.function.default(fun))

}



# @note
# \link[rpart]{rpart} is quite slow








#' @details
#' 
#' ## Dichotomize Multiple Predictors   
#' 
#' Function [m_rpartD] dichotomizes 
#' each predictor `X[,i]` based on the response \eqn{y} 
#' using function [rpartD].
#' Applying the multiple dichotomizing rules to a new set of predictors `newX`,
#' \itemize{
#' \item {A \link[base]{warning} message is produced, 
#' if at least one of the dichotomized predictors is all-`TRUE` or all-`FALSE`.}
#' \item {We do not check if more than one of the dichotomized predictors
#' are \link[base]{identical} to each other.
#' We take care of this situation in helper function [coef_dichotom]}
#' }
#' 
#' 
#' @returns 
#' 
#' ## Dichotomize Multiple Predictors
#' 
#' Function [m_rpartD] returns a \link[base]{function},
#' with a \link[base]{double} \link[base]{matrix} parameter `newX`. 
#' The argument for `newX` must have 
#' the same number of columns and the same column names as 
#' the input \link[base]{matrix} \eqn{X}.
#' The returned value of `m_rpartD(y,X)(newX)` is a 
#' \link[base]{logical} \link[base]{matrix}
#' with \link[base]{attributes}
#' \describe{
#' \item{`attr(,'cutoff')`}{
#' named \link[base]{double} \link[base]{vector}, 
#' the cutoff values for each predictor in `newX`}
#' }
#' 
#' @examples
#' ## Dichotomize Multiple Predictors
#' library(survival)
#' data(stagec, package = 'rpart') # see more details from ?rpart::stagec
#' nrow(stagec) # 146
#' (foo = with(stagec[1:100,], m_rpartD(y = Surv(pgtime, pgstat), X = cbind(age, g2, gleason))))
#' foo(as.matrix(stagec[-(1:100), c('age', 'g2', 'gleason')]))
#' @rdname rpartD
#' @export
m_rpartD <- function(y, X, check_degeneracy = TRUE, ...) {
  if (!is.matrix(X)) stop('`X` must be matrix')
  
  sq <- seq_len(dim(X)[2L]) # column sequence
  names(sq) <- dimnames(X)[[2L]]
  
  funs <- lapply(sq, FUN = function(i) rpartD(y = y, x = X[,i], check_degeneracy = FALSE, ...)) # a list of functions
  fbodies <- lapply(funs, FUN = body)
  
  calls <- lapply(sq, FUN = function(i) {
    cl <- fbodies[[i]]
    cl[[2L]] <- call('[', quote(newX), alist(a =)[[1L]], i)
    call('(', cl) # easier to read for human
  })
  
  cutoff <- vapply(fbodies, function(i) i[[3L]], FUN.VALUE = NA_real_)
  
  fun <- alist(newX = )
  fun[[2]] <- call(
    name = '{', 
    quote(if (!is.matrix(newX)) stop('input must be matrix')),
    call('if', call('!=', quote(dim(newX)[2L]), length(funs)), quote(stop())),
    call('if', call('!', call('identical', quote(dimnames(newX)[[2L]]), names(funs))), quote(stop())),
    # add: do not allow all-NA columns in `X`
    
    call('<-', quote(ret), as.call(c(list(name = quote(cbind)), calls))),
    
    if (check_degeneracy) quote(ms <- colMeans(ret, na.rm = TRUE)),
    if (check_degeneracy) quote(if (any(ms < .Machine$double.eps) || any(ms > 1 - .Machine$double.eps)) warning('Dichotomized value is all-0 or all-1')),

    # do not need to use `call('attr<-', ...)`
    call('<-', call('attr', quote(ret), which = 'cutoff'), cutoff),
    
    quote(return(ret))
  )
  
  return(as.function.default(fun))
  
}


