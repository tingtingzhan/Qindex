

#' @title Optimal Dichotomizing Predictors via Repeated Sample Splits
#' 
#' @description
#' To identify the optimal dichotomizing predictors using repeated sample splits.
#' 
#' @param formula,y,x \link[stats]{formula}, e.g., `y~X` or `y~x1+x2`.
#' Types of response \eqn{y} may be \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.
#' Candidate \link[base]{numeric} predictors \eqn{x}'s may be specified as the columns of one \link[base]{matrix} column, e.g., `y~X`; or as several \link[base]{vector} columns, e.g., `y~x1+x2`.
#' In helper functions, `x` is a \link[base]{numeric} \link[base]{vector}.
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param include (optional) \link[base]{language}, inclusion criteria. 
#' Default `(p1>.15 & p1<.85)` specifies a user-desired range of \eqn{p_1}
#' for the candidate dichotomizing predictors.
#' See explanation of \eqn{p_1} in section **Returns of Helper Functions**.
#' 
#' @param top positive \link[base]{integer} scalar, number of optimal dichotomizing predictors, default `1L`
#' 
#' @param nsplit,... additional parameters for function [rSplit]
#' 
#' @param id \link[base]{logical} \link[base]{vector} for helper function [split_dichotom], indices of training (`TRUE`) and test (`FALSE`) subjects 
#' @param ids (optional) \link[base]{list} of \link[base]{logical} \link[base]{vector}s for helper function [splits_dichotom], multiple copies of indices of repeated training-test sample splits.  
#' @param probs \link[base]{double} scalar for helper function [quantile.splits_dichotom], see \link[stats]{quantile}
#' 
#' 
#' @details 
#' 
#' Function [optimSplit_dichotom] selects the optimal dichotomizing predictors via repeated sample splits. Specifically, 
#' 
#' \enumerate{
#' \item Generate multiple, i.e., repeated, training-test sample splits (via [rSplit])
#' \item For each candidate predictor \eqn{x_i}, find the ***median-split-dichotomized regression model*** based on the repeated sample splits, see details in section **Details on Helper Functions**
#' \item Limit the selection of the candidate predictors \eqn{x}'s to a user-desired range of \eqn{p_1} of the split-dichotomized regression models, see explanations of \eqn{p_1} in section **Returns of Helper Functions**
#' \item Rank the candidate predictors \eqn{x}'s by the decreasing order of the \link[base]{abs}olute values of the regression coefficient estimate of the median-split-dichotomized regression models.  On the top of this rank are the ***optimal dichotomizing predictors***.
#' }
#' 
#' @returns 
#' Function [optimSplit_dichotom] returns an object of \link[base]{class} `'optimSplit_dichotom'`, which is a \link[base]{list} of dichotomizing \link[base]{function}s, 
#' with the input `formula` and `data` as additional \link[base]{attributes}.
#' 
#' @examples 
#' # see ?`Qindex-package`
#' @rdname optimSplit_dichotom
#' @export
optimSplit_dichotom <- function(
    formula, data,
    include = quote(p1 > .15 & p1 < .85), # ?devtools::check warning
    top = 1L,
    nsplit,
    ...
) {
  
  cl <- match.call()
  y <- eval(formula[[2L]], envir = data)
  ids <- rSplit(y, nsplit = nsplit, ...) # using same split for all predictors
  
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
    tmp <- splits_dichotom(y = y, x = X[,p], ids = ids)
    quantile.splits_dichotom(tmp, probs = .5)
  })
  mssd <- do.call(what = Map, args = c(list(f = c), lapply(tmp, FUN = function(i) attributes(i)[c('rule', 'text', 'p1', 'coef')])))
  
  ret_rule <- mapply(FUN = function(rule, nm) {
    # rule = mssd$rule[[1L]]; nm = colnames(X)[[1L]]
    r0 <- as.list.function(rule)
    r0[[1L]] <- if (is.symbol(formula[[3L]])) {
      call(name = '[', formula[[3L]], alist(i =)[[1L]], nm)
    } else stop('should be even easier')
    return(as.function.default(r0))
  }, rule = mssd$rule, nm = colnames(X))
  
  names(ret_rule) <- colnames(X)
  
  if (!missing(include)) include <- substitute(include)
  id_include <- eval(call(name = 'with.default', data = quote(mssd), expr = include))
  mssd$coef[!id_include] <- NA_real_
  
  id_top <- order(abs(mssd$coef), decreasing = TRUE)[seq_len(top)]
  if (anyNA(mssd$coef[id_top])) stop('Decrease `top` (containing coef\'s which do not satisfy `include`)')

  ret <- ret_rule[id_top]
  attr(ret, which = 'formula') <- formula
  attr(ret, which = 'data') <- data
  class(ret) <- 'optimSplit_dichotom'
  return(ret)
  
}


#' @export
`[.optimSplit_dichotom` <- function(x, i) {
  ret <- unclass(x)[i] # NOT `[[`
  attr(ret, which = 'formula') <- attr(x, which = 'formula', exact = TRUE)
  attr(ret, which = 'data') <- attr(x, which = 'data', exact = TRUE)
  class(ret) <- 'optimSplit_dichotom'
  return(ret)
}




#' @export
print.optimSplit_dichotom <- function(x, ...) {
  lapply(x, FUN = function(rule) { # (rule = x[[1L]])
    out <- body(rule)[[2L]][[3L]][[2L]]
    out[[2L]] <- as.list.function(rule)[[1L]]
    print(out)
    return(invisible())
  })
  return(invisible())
}



#' @title Regression Models with Optimal Dichotomizing Predictors
#' 
#' @description
#' Regression models with optimal dichotomizing predictor(s), used either as boolean or continuous predictor(s).
#' 
#' @param object an [optimSplit_dichotom] object
#' @param formula (optional) \link[stats]{formula} to specify the response in test data. If missing, the model formula of training data is used
#' @param newdata (optional) test \link[base]{data.frame}, candidate \link[base]{numeric} predictors \eqn{x}'s must have the same \link[base]{name} and \link[base]{dim}ension as the training data. If missing, the training data is used
#' @param boolean \link[base]{logical} scalar, whether to use the *dichotomized* predictor (default, `TRUE`), or the continuous predictor (`FALSE`)
#' @param ... additional paramters, currently not in use
#' 
#' @returns
#' Function [predict.optimSplit_dichotom] returns a \link[base]{list} of regression models, \link[survival]{coxph} model for \link[survival]{Surv} response, \link[stats]{glm} for \link[base]{logical} response, and \link[stats]{lm} model for \link[base]{numeric} response.
#' 
#' @examples
#' # see ?`Qindex-package`
#' @export predict.optimSplit_dichotom
#' @export
predict.optimSplit_dichotom <- function(
    object, 
    formula = attr(object, which = 'formula', exact = TRUE),
    newdata = attr(object, which = 'data', exact = TRUE),
    boolean = TRUE,
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
    
    if (boolean) {
      bool_ <- with(newdata, do.call(rule, args = list(nm0)))
      # why do I have to do this?? why neither of following does not work?
      #eval(rule(), envir = newdata)
      #with(newdata, rule())
      nm_ <- paste0(deparse1(nm0), attr(bool_, which = 'text', exact = TRUE))
      assign(x = nm_, value = bool_)
    } else { # using continuous (i.e., not dichotomized) predictor
      nm_ <- deparse1(nm0)
      assign(x = nm_, value = eval(nm0, envir = newdata))
    }
    
    fom_ <- eval(call('~', quote(y), as.symbol(nm_)))
    
    suppressWarnings(if (inherits(y, what = 'Surv')) {
      do.call('coxph', args = list(formula = fom_))
    } else if (is.logical(y) || all(y %in% c(0, 1))) {
      do.call('glm', args = list(formula = fom_, family = quote(binomial(link = 'logit'))))
    } else do.call('lm', args = list(formula = fom_)))
  })
  
  return(ret)
  
}















#' @section Details on Helper Functions:
#' 
#' ## Split-Dichotomized Regression Model
#' 
#' Helper function [split_dichotom] performs a univariable regression model on the test set with a dichotomized predictor, using a dichotomizing rule determined by a recursive partitioning of the training set. 
#' Specifically, given a training-test sample split,
#' \enumerate{
#' \item find the *dichotomizing rule* \eqn{\mathcal{D}} of the predictor \eqn{x_0} given the response \eqn{y_0} in the training set (via [rpartD]);
#' \item fit a univariable regression model of the response \eqn{y_1} with the dichotomized predictor \eqn{\mathcal{D}(x_1)} in the test set.
#' }
#' Currently the Cox proportional hazards (\link[survival]{coxph}) regression for \link[survival]{Surv} response, logistic (\link[stats]{glm}) regression for \link[base]{logical} response and linear (\link[stats]{lm}) regression for \link[stats]{gaussian} response are supported.
#' 
#' @section Returns of Helper Functions: 
#' 
#' Helper function [split_dichotom] returns a split-dichotomized regression model, which is either a Cox proportional hazards (\link[survival]{coxph}), a logistic (\link[stats]{glm}), or a linear (\link[stats]{lm}) regression model, with additional \link[base]{attributes}
#' 
#' \describe{
#' \item{`attr(,'rule')`}{\link[base]{function}, dichotomizing rule \eqn{\mathcal{D}} based on the training set}
#' \item{`attr(,'text')`}{\link[base]{character} scalar, human-friendly description of \eqn{\mathcal{D}}}
#' \item{`attr(,'p1')`}{\link[base]{double} scalar, \eqn{p_1 = \text{Pr}(\mathcal{D}(x_1)=1)}}
#' \item{`attr(,'coef')`}{\link[base]{double} scalar, univariable regression coefficient estimate of \eqn{y_1\sim\mathcal{D}(x_1)}}
#' }
#' 
#' @importFrom survival coxph
#' @importFrom stats lm glm binomial
#' @rdname optimSplit_dichotom
#' @export
split_dichotom <- function(y, x, id, ...) {
  
  # id: training set
  # !id: test set
  rule <- rpartD(y = y[id], x = x[id], check_degeneracy = TRUE)
  dtest <- tryCatch(data.frame(y = y[!id], dx = rule(x[!id])), warning = identity)
  
  if (inherits(dtest, what = 'warning')) {
    # exception
    mtest <- logical()
    attr(mtest, which = 'coef') <- NA_real_
    return(mtest)
  }
  
  suppressWarnings(mtest <- if (inherits(y, what = 'Surv')) {
    coxph(formula = y ~ dx, data = dtest)
  } else if (is.logical(y) || all(y %in% c(0, 1))) {
    glm(formula = y ~ dx, family = binomial(link = 'logit'), data = dtest)
  } else lm(formula = y ~ dx, data = dtest))
  
  cf_test <- mtest$coefficients[length(mtest$coefficients)]
  attr(mtest, which = 'rule') <- rule
  attr(mtest, which = 'text') <- attr(dtest$dx, which = 'text', exact = TRUE)
  attr(mtest, which = 'p1') <- mean.default(dtest$dx, na.rm = TRUE)
  attr(mtest, which = 'coef') <- if (is.finite(cf_test)) unname(cf_test) else NA_real_
  
  class(mtest) <- c('split_dichotom', class(mtest))
  return(mtest)
  
}


#' @section Details on Helper Functions:
#' 
#' ## Split-Dichotomized Regression Models based on Repeated Training-Test Sample Splits
#' 
#' Helper function [splits_dichotom] fits multiple split-dichotomized regression models [split_dichotom] on the response \eqn{y} and predictor \eqn{x}, based on each copy of the repeated training-test sample splits.
#' 
#' @section Returns of Helper Functions: 
#' 
#' Helper function [splits_dichotom] returns a \link[base]{list} of split-dichotomized regression models ([split_dichotom]).
#' 
#' @rdname optimSplit_dichotom
#' @export
splits_dichotom <- function(y, x, ids = rSplit(y, ...), ...) {
  ret <- lapply(ids, FUN = function(id) split_dichotom(y, x, id = id, ...))
  class(ret) <- c('splits_dichotom', class(ret))
  return(ret)
}
  




#' @section Details on Helper Functions:
#' 
#' ## Quantile of Split-Dichotomized Regression Models
#' 
#' Helper function [quantile.splits_dichotom] is a method dispatch of the S3 generic function \link[stats]{quantile} on [splits_dichotom] object.
# finds the \link[stats]{quantile}
# of the univariable regression coefficient (i.e., effect size) of a dichotomized predictor,
# based on repeated given training-test sample splits.
#' Specifically,
#' 
#' \enumerate{
#' \item {collect the univariable regression coefficient estimate from each one of the split-dichotomized regression models;}
#' \item {find the nearest-even (i.e., `type = 3`) \link[stats]{quantile} of the coefficients from Step 1. By default, we use the \link[stats]{median} (i.e., `prob = .5`);}
#' \item {the split-dichotomized regression model corresponding to the selected coefficient quantile in Step 2, is returned.}
#' }
#' 
#' @section Returns of Helper Functions: 
#' Helper function [quantile.splits_dichotom] returns a split-dichotomized regression model ([split_dichotom]).
#' 
#' @importFrom stats quantile
#' @rdname optimSplit_dichotom
#' @export quantile.splits_dichotom
#' @export
# old name [quantile_split_dichotom]
quantile.splits_dichotom <- function(x, probs = .5, ...) {
  cf <- vapply(x, FUN = attr, which = 'coef', exact = TRUE, FUN.VALUE = NA_real_)
  medianID <- which(cf == quantile(cf, probs = probs, type = 3L, na.rm = TRUE))[1L]
  return(x[[medianID]])
}







