
#' @title Qindex via Cross Validation
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param k \link[base]{integer} scalar
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @examples
#' data(wrobel_lung, package = 'spatstat.grouped.data')
#' library(spatstat.grouped)
#' lungQp = grouped_ppp(hladr ~ OS | patient_id/image_id, f_sum_ = 'min', data = wrobel_lung) |> 
#'  aggregate_quantile(by = ~ patient_id, probs = seq.int(from = .05, to = .95, by = .01))
#' 
#' lung_CV_QI = Qindex_cv(OS ~ hladr.quantile, data = lungQp, k = 10L, nonlinear = FALSE)
#' head(lung_CV_QI)   
#' boxplot(QI_notAdj ~ folds., data = lung_CV_QI) 
#' boxplot(QI_sgnAdj ~ folds., data = lung_CV_QI) 
#' boxplot(QI_globalAdj ~ folds., data = lung_CV_QI) 
#' table(attr(lung_CV_QI, 'sign'))
#' library(survival)
#' summary(coxph(OS ~ QI_notAdj, data = lung_CV_QI))
#' summary(coxph(OS ~ QI_sgnAdj, data = lung_CV_QI))
#' summary(coxph(OS ~ QI_globalAdj, data = lung_CV_QI))
#' 
#' @importFrom caret createFolds
#' @importFrom gam.matrix gam_matrix predict.gam_matrix cor_xy.gam_matrix
#' @export
Qindex_cv <- function(formula, data, k, ...) { 

  fld <- createFolds(y = seq_len(.row_names_info(data, type = 2L)), k = k, list = TRUE, returnTrain = FALSE)
  
  dataout <- data
  dataout$folds. <- NA_integer_
  dataout$QI_notAdj <- dataout$QI_sgnAdj <- NA_real_
  sgn <- integer(length = k)
  
  for (i in seq_along(fld)) {
    id <- fld[[i]]
    dataout$folds.[id] <- i
    d0 <- data[-id, , drop = FALSE] # training set
    d1 <- data[id, , drop = FALSE] # test set
    m <- gam_matrix(formula = formula, data = d0, ...)
    dataout$QI_sgnAdj[id] <- predict.gam_matrix(m, newdata = d1, sign_adjusted = TRUE)
    dataout$QI_notAdj[id] <- predict.gam_matrix(m, newdata = d1, sign_adjusted = FALSE)
    sgn[i] <- cor_xy.gam_matrix(m) |> sign()
  }
  
  m_global <- gam_matrix(formula = formula, data = data, ...)
  sign_global <- cor_xy.gam_matrix(m_global) |> sign()
  dataout$QI_globalAdj <- dataout$QI_notAdj * sign_global
  
  attr(dataout, which = 'sign') <- sgn
  
  return(dataout)
  
}
  


