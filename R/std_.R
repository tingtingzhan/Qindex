
#' @title Alternative Standardization Methods
#' 
#' @description
#' Alternative standardize using \link[stats]{median}, \link[stats]{IQR} and \link[stats]{mad}.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}
#' 
#' @param na.rm \link[base]{logical} scalar, 
#' see \link[stats]{quantile}, \link[stats]{median} and \link[stats]{mad}.
#' Default `TRUE`
#' 
#' @param ... additional parameters of \link[stats]{quantile} and/or \link[stats]{mad}
#' 
#' @return 
#' 
#' ## Standardize using \link[stats]{median} and \link[stats]{IQR}
#' Function [std_IQR] returns a \link[base]{numeric} \link[base]{vector} of the same length as `x`.
#' 
#' ## Standardize using \link[stats]{median} and \link[stats]{mad}
#' Function [std_mad] returns a \link[base]{numeric} \link[base]{vector} of the same length as `x`.
#' 
#' @examples
#' std_IQR(rnorm(20))
#' std_mad(rnorm(20))
#' 
#' @keywords internal
#' @importFrom stats quantile
#' @name std_
#' @export
std_IQR <- function(x, na.rm = TRUE, ...) {
  qs <- quantile(x, probs = c(.25, .5, .75), na.rm = na.rm, ...)
  (x - qs[2L]) / (qs[3L] - qs[1L])
}



#' @importFrom stats median.default mad
#' @rdname std_
#' @export
std_mad <- function(x, na.rm = TRUE, ...) {
  m <- median.default(x, na.rm = na.rm)
  mad_ <- mad(x, center = m, na.rm = na.rm, ...)
  (x - m) / mad_
}




