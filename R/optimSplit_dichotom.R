

#' @title Optimal Dichotomizing Predictors via Repeated Sample Splits
#' 
#' @description
#' 
#' Functions explained in this documentation are,
#' 
#' \describe{
#' 
#' \item{[optimSplit_dichotom]}{
#' to identify the optimal dichotomizing predictors using repeated sample splits.}
#' 
#' \item{[split_dichotom]}{
#' a helper function to perform a univariable regression model on the test set 
#' with a dichotomized predictor,
#' using a dichotomizing rule determined 
#' by a recursive partitioning of the training set.}
#' 
#' \item{[quantile_split_dichotom]}{
#' a helper function to locate a quantile of multiple [split_dichotom] objects, 
#' based on the estimated univariable regression coefficient.}
#' 
#' }
#' 
#' 
#' 
#' 
#' @param formula \link[stats]{formula}. 
#' Left-hand-side is the \link[base]{name} of 
#' a \link[survival]{Surv}, \link[base]{logical}, or \link[base]{double} response \eqn{y}.
#' Right-hand-side is the candidate \link[base]{numeric} predictors in `data`, 
#' given either as the \link[base]{name} of a \link[base]{numeric} \link[base]{matrix} column
#' (e.g., `y ~ X`), 
#' or as the names of several \link[base]{numeric} \link[base]{vector} columns
#' (e.g., `y ~ x1 + x2 + x3`)
#' 
#' @param data \link[base]{data.frame}, containing the response and predictors in `formula`
#' 
#' @param include \link[base]{language} object, 
#' inclusion criteria for the optimal dichotomizing predictors. 
#' A suggested choice is `(highX>.15 & highX<.85)`
#' to guarantee a user-desired range of proportions in `highX`.
#' See explanation of `highX` in helper function [split_dichotom].
#' 
#' @param top positive \link[base]{integer} scalar, number of optimal dichotomizing predictors, default `1L`
#' 
#' @param nsplit,... additional parameters for function [rSplit]
#' 
#' @param y (for helper functions) 
#' a \link[survival]{Surv} object, a \link[base]{logical} \link[base]{vector}, 
#' or a \link[base]{double} \link[base]{vector}, the response \eqn{y}
#' 
#' @param x (for helper functions) 
#' \link[base]{numeric} \link[base]{vector}, a single predictor \eqn{x}
#' 
#' @param index (for helper function [split_dichotom]) 
#' \link[base]{logical} \link[base]{vector},
#' indices of training and test set. 
#' `TRUE` elements indicate training subjects and 
#' `FALSE` elements indicate test subjects.
#' 
#' @param indices (optional, for helper function [quantile_split_dichotom]) 
#' a \link[base]{list} of \link[base]{logical} \link[base]{vector}s,
#' the indices of multiple training-test sample splits.  
#' Default value is provided by function [rSplit].
#' 
#' @param probs (for helper function [quantile_split_dichotom]) 
#' \link[base]{double} scalar, see \link[stats]{quantile}
#' 
#' 
#' @details 
#' 
#' Function [optimSplit_dichotom] selects the optimal dichotomizing predictors via repeated sample splits.
#' Specifically,
#' 
#' \enumerate{
#' 
#' \item Generate multiple training-test sample splits using function [rSplit]
#' 
#' \item For each candidate predictor, 
#' find the median [split_dichotom] 
#' (using helper function [quantile_split_dichotom]) 
#' of the multiple sample splits from Step 1.
#' 
#' \item (Optional) limit the selection in a subset of the candidate predictors.
#' Typically, we would prefer to guarantee
#' a user-desired range of `highX` 
#' (see explanations on `highX` in section **Returns of Helper Functions**).
#' A suggested choice is `(highX>.15 & highX<.85)`.
#' 
#' \item Rank the candidate predictors, from either Step 2 or Step 3, 
#' by the decreasing order of the \link[base]{abs}olute values of 
#' the estimated univariable regression coefficients of the corresponding [split_dichotom] objects.
#' 
#' }
#' 
#' The *optimal dichotomizing predictors* are the ones
#' with the largest \link[base]{abs}olute values of 
#' the estimated univariable regression coefficients 
#' of the corresponding [split_dichotom] objects.
#' 
#' 
#' 
#' @returns 
#' Function [optimSplit_dichotom] returns an object of \link[base]{class} `'optimSplit_dichotom'`,
#' which is a \link[base]{list} of dichotomizing \link[base]{function}s
#' of \link[base]{length}-`top`.
# which contains the response, 
# and only the optimal dichotomizing predictors out of all candidate predictors.
# Other variables in `data`, which are not specified in `formula`, are retained.
# In addition, the dichotomized values of the optimal dichotomizing predictors,
# according to their respective dichotomizing rules, are also included.
#' The returned value has \link[base]{attributes},
#' \describe{
#' \item{`attr(,'formula')`}{
#' input \link[stats]{formula}}
#' \item{`attr(,'data')`}{
#' input \link[base]{data.frame}}
# \item{`attr(,'id_top')`}{
# positive \link[base]{integer} scalar or \link[base]{vector},
# the indices of the optimal dichotomizing predictors out of all candidate predictors.}
# \item{`attr(,'top')`}{
# a diagnostic \link[base]{data.frame} of 
# the median [split_dichotom]s of each of the optimal dichotomizing predictors,
# with columns 
# \describe{
# \item{`$cutoff`}{the cutoff threshold, identified in the training set}
# \item{`$highX`}{
# proportion of the dichotomizing predictors 
# greater-than or greater-than-or-equal-to the cutoff threshold, in the test set}
# \item{`$coef`}{
# the estimated univariable regression coefficient of 
# the dichotomized predictor, in the test set}
# }
# }
#' }
#' 
#' 
#' @rdname optimSplit_dichotom
#' @export
optimSplit_dichotom <- function(
    formula, data,
    #include = (highX > .15 & highX < .85), # ?devtools::check warning
    include,
    top = 1L,
    nsplit,
    ...
) {
  
  cl <- match.call()
  y <- eval(formula[[2L]], envir = data)
  indices <- rSplit(y, nsplit = nsplit, ...) # using same split for all predictors
  
  if (is.symbol(formula[[3L]])) {
    X <- eval(formula[[3L]], envir = data) # 'matrix' of predictors
  } else {
    fom <- eval(call(name = '~', call(name = '+', formula[[3L]], quote((-1)))))
    X <- as.matrix.data.frame(model.frame.default(formula = fom, data = data, na.action = na.pass))
    # ?stats::model.matrix.default does not have parameter `na.action`
  }
  
  if (!is.numeric(X) || !is.matrix(X)) stop('predictors to be dichotomized must be numeric')
  #if (anyNA(X)) # it's okay now!
  
  tmp <- lapply(seq_len(dim(X)[2L]), FUN = function(p) {
    quantile_split_dichotom(y = y, x = X[,p], indices = indices, probs = .5)
  })
  mssd <- do.call(what = Map, args = c(list(f = c), lapply(tmp, FUN = function(i) attributes(i)[c('rule', 'cutoff', 'text', 'highX', 'coef')])))
  
  ret_rule <- mapply(FUN = function(rule, nm) {
    # rule = mssd$rule[[1L]]; nm = colnames(X)[[1L]]
    r0 <- as.list.function(rule)
    r0[[1L]] <- if (is.symbol(formula[[3L]])) {
      call(name = '[', formula[[3L]], alist(i =)[[1L]], nm)
    } else stop('should be even easier')
    return(as.function.default(r0))
  }, rule = mssd$rule, nm = colnames(X))
  
  if (!missing(include)) {
    id_excl <- !eval(call(name = 'with.default', data = quote(mssd), expr = substitute(include)))
    mssd$coef[id_excl] <- NA_real_
  }
  
  id_top <- order(abs(mssd$coef), decreasing = TRUE)[seq_len(top)]
  if (anyNA(mssd$coef[id_top])) stop('Decrease `top` (containing coef\'s which do not satisfy `include`)')

  ret <- ret_rule[id_top]
  
  # to be removed!!!  
  #nm_top <- colnames(X)[id_top]
  #names(ret) <- nm_top
  # end of to be removed!!!
  
  attr(ret, which = 'formula') <- formula
  attr(ret, which = 'data') <- data
  #attr(ret, which = 'call') <- cl
  class(ret) <- 'optimSplit_dichotom'
  return(ret)
  
}

# @export
#print.optimSplit_dichotom <- function(x, ...) {
#  .Defunct(msg = 'no need for now')
#  mapply(FUN = function(rule, nm) {
#    ret <- body(rule)[[2L]][[3L]][[2L]]
#    ret[[2L]] <- as.symbol(nm)
#    print(ret)
#    return(invisible())
#  }, rule = x, nm = names(x))
#  return(invisible())
#}



#' @title predict.optimSplit_dichotom
#' 
#' @description
#' ..
#' 
#' @param object an [optimSplit_dichotom] object
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param newdata (optional) test \link[base]{data.frame}.
#' If missing, the training data is used
#' 
#' @param ... ..
#' 
#' @returns
#' Function [predict.optimSplit_dichotom] returns ..
#' 
#' @export predict.optimSplit_dichotom
#' @export
predict.optimSplit_dichotom <- function(
    object, 
    formula = attr(object, which = 'formula', exact = TRUE),
    newdata = attr(object, which = 'data', exact = TRUE), 
    ...
) {
  
  y <- eval(formula[[2L]], envir = newdata)
  
  # KEEP FOR NOW
  #nm <- if (is.symbol(formula[[3L]])) {# 'matrix' X
  #  sprintf(fmt = '%s[,\'%s\']', as.character(formul[[3L]]), names(object))
  #} else names(object)
  # I need to make sure `!is.symbol(formula[[3L]])` works fine too
  
  ret <- lapply(object, FUN = function(rule) {
    # rule = object[[1L]]
    
    nm0 <- as.list.function(rule)[[1L]]
    
    bool_ <- with(newdata, do.call(rule, args = list(nm0)))
    # why do I have to do this?? why neither of following does not work?
    #eval(rule(), envir = newdata)
    #with(newdata, rule())
    
    nm_ <- paste0(deparse1(nm0), attr(bool_, which = 'text'))
    assign(x = nm_, value = bool_)
    fom_ <- eval(call('~', quote(y), as.symbol(nm_)))
    if (inherits(y, what = 'Surv')) {
      suppressWarnings(do.call('coxph', args = list(formula = fom_)))
    } else if (is.logical(y) || all(y %in% c(0, 1))) {
      suppressWarnings(do.call('glm', args = list(formula = fom_, family = binomial(link = 'logit'))))
    } else if (is.vector(y, mode = 'numeric')) {
      suppressWarnings(do.call('lm', args = list(formula = fom_)))
    }
  })
  
  return(ret)
  
}















#' @section Details on Helper Functions:
#' 
#' ## Univariable regression model with a dichotomized predictor
#' 
#' Helper function [split_dichotom] performs a univariable regression model on the test set 
#' with a dichotomized predictor,
#' using a dichotomizing rule determined 
#' by a recursive partitioning of the training set. 
#' Currently the Cox proportional hazards (\link[survival]{coxph}) regression for \link[survival]{Surv} response, 
#' logistic (\link[stats]{glm}) regression for \link[base]{logical} response and 
#' linear (\link[stats]{lm}) regression for \link[stats]{gaussian} response
#' are supported.
#' Specifically, given a training-test sample split,
#' 
#' \enumerate{
#' \item find the dichotomizing rule of the response \eqn{y} 
#' given the predictor \eqn{x}, using function [rpartD], in the training set
#' \item dichotomize the predictor \eqn{x} using the rule identified in Step 1, 
#' in the test set.
#' \item run a univariable regression model on the response \eqn{y} 
#' on the dichotomized predictor from Step 2, in the test set.
#' }
#' 
#' 
#' 
#' 
#' @section Returns of Helper Functions: 
#' 
#' Helper function [split_dichotom], as well as helper function [quantile_split_dichotom], returns 
#' a Cox proportional hazards (\link[survival]{coxph}), 
#' or a logistic (\link[stats]{glm}), 
#' or a linear (\link[stats]{lm}) 
#' regression model, 
#' with additional \link[base]{attributes}
#' 
#' \describe{
#' \item{`attr(,'rule')`}{\link[base]{function}, 
#' the dichotomizing rule based on the training set}
#' \item{`attr(,'cutoff')`}{\link[base]{numeric} scalar, 
#' the cutoff threshold based on the training set}
#' \item{`attr(,'highX')`}{\link[base]{double} scalar, 
#' proportion of \link[base]{numeric} predictor \eqn{x}, in the test set, which is greater-than or greater-than-or-equal-to
#' the cutoff threshold `attr(, 'cutoff')`}
#' \item{`attr(,'coef')`}{\link[base]{double} scalar, 
#' the estimated univariable regression coefficient of the dichotomized predictor in the test set}
#' }
#' 
#' 
#' 
#' @importFrom survival coxph
#' @importFrom stats lm glm binomial
#' @rdname optimSplit_dichotom
#' @export
split_dichotom <- function(y, x, index, ...) {
  
  # index: training set
  # !index: test set
  branch <- rpartD(y = y[index], x = x[index], check_degeneracy = TRUE)
  dtest <- tryCatch(data.frame(y = y[!index], high = branch(x[!index])), warning = identity)
  
  if (inherits(dtest, what = 'warning')) {
    # exception
    mtest <- logical()
    attr(mtest, which = 'coef') <- NA_real_
    return(mtest)
  }
  
  mtest <- if (inherits(y, what = 'Surv')) {
    suppressWarnings(coxph(formula = y ~ high, data = dtest))
  } else if (is.logical(y) || all(y %in% c(0, 1))) {
    suppressWarnings(glm(formula = y ~ high, family = binomial(link = 'logit'), data = dtest))
  } else if (is.vector(y, mode = 'numeric')) {
    suppressWarnings(lm(formula = y ~ high, data = dtest))
  }
  
  cf_test <- mtest$coefficients[length(mtest$coefficients)]
  attr(mtest, which = 'coef') <- if (is.finite(cf_test)) unname(cf_test) else NA_real_
  attr(mtest, which = 'highX') <- mean.default(dtest$high, na.rm = TRUE)
  attr(mtest, which = 'rule') <- branch
  attr(mtest, which = 'cutoff') <- attr(dtest$high, which = 'cutoff', exact = TRUE)
  attr(mtest, which = 'text') <- attr(dtest$high, which = 'text', exact = TRUE)
  
  # class(mtest) <- c('split_dichotom', class(mtest))
  return(mtest)
  
}











#' @section Details on Helper Functions:
#' 
#' ## Quantile of [split_dichotom] objects
#' 
#' Helper function [quantile_split_dichotom] finds the \link[stats]{quantile}
#' of the univariable regression coefficient (i.e., effect size) of a dichotomized predictor,
#' based on multiple given training-test sample splits.
#' Specifically,
#' 
#' \enumerate{
#' \item {for each training-test sample split, 
#' fit the univariable regression model based on the dichotomized predictor, 
#' using helper function [split_dichotom]}
#' \item {finds the nearest-even (`type = 3`) \link[stats]{quantile}
#' of the estimated univariable regression coefficients obtained in Step 1, 
#' based on the user-specified probability `prob`}
#' }
#' 
#' The [split_dichotom] object from Step 1, 
#' whose estimated univariable regression coefficient equals to 
#' the specified quantile identified in Step 2,
#' is referred to as the quantile of [split_dichotom] objects 
#' based on the multiple given training-test sample splits.
#' 
#' @importFrom stats quantile
#' @rdname optimSplit_dichotom
#' @export
quantile_split_dichotom <- function(y, x, indices = rSplit(y, ...), probs = .5, ...) {
  
  tmp <- lapply(indices, FUN = function(index) split_dichotom(y, x, index = index, ...))
  
  cf <- vapply(tmp, FUN = attr, which = 'coef', exact = TRUE, FUN.VALUE = NA_real_)
  
  medianID <- which(cf == quantile(cf, probs = probs, type = 3L, na.rm = TRUE))[1L]
  return(tmp[[medianID]])
  
}







