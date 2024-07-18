
#' @title Bootstrap Indices
#' 
#' @description
#' Generate a series of \link[boot]{boot}strap indices.
#' 
#' @param n positive \link[base]{integer} scalar, sample size \eqn{n}
#' 
#' @param R positive \link[base]{integer} scalar, number of bootstrap replicates \eqn{R}
#' 
#' @returns 
#' Function [bootid] returns a length-\eqn{R} \link[base]{list} of 
#' positive \link[base]{integer} \link[base]{vector}s.
#' Each element is the length-\eqn{n} indices of each bootstrap sample.
#' 
#' @details
#' Function [bootid] generates the same bootstrap indices as 
#' those generated from the default options of function \link[boot]{boot} 
#' (i.e., `sim = 'ordinary'` and `m = 0`), 
#' given the same \link[base]{Random} seed.  
#' 
#' See details in `boot:::index.array` and `boot:::ordinary.array`.
#' 
#' @examples
#' set.seed(1345); boot::boot(data = 1:10, statistic = function(data, ind) ind, R = 3L)[['t']]
#' set.seed(1345); bootid(10L, R = 3L) # same copies of indices
#' 
#' @keywords internal
#' @export
bootid <- function(n, R) {
  ret <- sample.int(n = n, size = n * R, replace = TRUE)
  dim(ret) <- c(R, n)
  lapply(seq_len(R), FUN = function(i) ret[i,])
}

