
#' @title Stratified Random Split Sampling
#' 
#' @description 
#' Random split sampling, stratified based on the type of the response.
#' 
#' @param y a \link[base]{double} \link[base]{vector}, 
#' a \link[base]{logical} \link[base]{vector},
#' a \link[base]{factor}, 
#' or a \link[survival]{Surv} object, 
#' response \eqn{y}
#' 
#' @param stratify \link[base]{logical} scalar, 
#' whether stratification based on response \eqn{y} needs to be implemented, default `TRUE`
#' 
#' @param nsplit positive \link[base]{integer} scalar, number of \link[base]{replicate}s of random splits to be performed
#' 
#' @param s_ratio \link[base]{double} scalar between 0 and 1, 
#' split ratio, i.e., percentage of training subjects \eqn{p}, default `.8`
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' Function [rSplit] performs random split sampling, 
#' with or without stratification. Specifically,
#' 
#' \itemize{
#' 
#' \item If `stratify = FALSE`, 
#' or if we have a \link[base]{double} response \eqn{y},
#' then split the sample into a training and a test set by odds \eqn{p/(1-p)}, without stratification.
#' 
#' \item Otherwise, split a \link[survival]{Surv} response \eqn{y}, stratified by its censoring status.
#' Specifically, 
#' split subjects with observed event into a training and a test set by odds \eqn{p/(1-p)},
#' and split the censored subjects into a training and a test set by odds \eqn{p/(1-p)}.
#' Then combine the training sets from subjects with observed events and censored subjects,
#' and combine the test sets from subjects with observed events and censored subjects.
#' 
#' \item Otherwise, split a \link[base]{logical} response \eqn{y}, stratified by itself.
#' Specifically, 
#' split the subjects with `TRUE` response into a training and a test set by odds \eqn{p/(1-p)},
#' and split the subjects with `FALSE` response into a training and a test set by odds \eqn{p/(1-p)}.
#' Then combine the training sets, and the test sets, in a similar fashion as described above.
#' 
#' \item Otherwise, split a \link[base]{factor} response \eqn{y}, stratified by its \link[base]{levels}.
#' Specifically, 
#' split the subjects in each level of \eqn{y} into a training and a test set by odds \eqn{p/(1-p)}.
#' Then combine the training sets, and the test sets, from all levels of \eqn{y}.
#' 
#' }
#' 
#' 
#' @returns 
#' Function [rSplit] returns a length-`nsplit` \link[base]{list} of 
#' \link[base]{logical} \link[base]{vector}s.
#' In each \link[base]{logical} \link[base]{vector}, 
#' the `TRUE` elements indicate training subjects and 
#' the `FALSE` elements indicate test subjects.
#' 
#' @note
#' `caTools::sample.split` is not what we need.
#' 
#' @examples
#' rSplit(y = rep(c(TRUE, FALSE), times = c(20, 30)), nsplit = 3L)
#' 
#' @seealso \link[base]{split}, `caret::createDataPartition`
#' @keywords internal
#' @export 
rSplit <- function(y, nsplit, stratify = TRUE, s_ratio = .8, ...) {
  
  if (anyNA(y)) stop('do not allow missingness in the response, for now')
  
  n <- length(y) 
  ret0 <- rep(FALSE, times = n)
  # works correctly for ?survival::Surv object, via ?survival:::length.Surv
  
  if (!stratify) {
    # no stratification
    idx <- list(seq_len(n))
    
  } else if (inherits(y, what = 'Surv')) {
    # stratify by censoring status
    if (dim(y)[2L] == 3L) stop('3-col Surv response not supported yet')
    xevent <- as.logical(y[,2L])
    idx <- list(
      which(!xevent), # 'integer' indices of censored events
      which(xevent) # 'integer' indices of observed events
    )

  } else if (is.logical(y) || all(y %in% c(0, 1))) {
    # stratify by the binary response
    y <- as.logical(y)
    idx <- list(
      which(!y), # 'integer' indices of non-responder
      which(y) # 'integer' indices of responder
    )

  } else if (is.factor(y)) {
    # stratify by the levels
    idx <- lapply(seq_along(attr(y, which = 'levels', exact = TRUE)), FUN = function(i) {
      which(unclass(y) == i)
    })
    
  } else if (is.vector(y, mode = 'numeric')) { # must after `binary` if
    # no stratification
    idx <- list(seq_len(n))
    
  } else stop('unsupported class: ', class(y)[1L])
  
  replicate(n = nsplit, expr = {
    idx_train <- lapply(idx, FUN = function(id) {
      sample(id, size = floor(length(id) * s_ratio), replace = FALSE)
    })
    #train <- sort.int(unlist(idx_train, use.names = FALSE))
    #list(train = train, test = sort.int(setdiff(seq_len(n), y = train)))
    train <- unlist(idx_train, use.names = FALSE)
    ret <- ret0
    ret[train] <- TRUE
    ret
  }, simplify = FALSE)
  
}

