

#' @title Continuous and Dichotomized Index Predictors Based on Distribution Quantiles
#' 
#' @description
#' Continuous and dichotomized index predictors based on distribution quantiles.
#' 
#' @section Step 1. Cluster-Specific Sample Quantiles [clusterQp]:
#'
#' Function [clusterQp] calculates user-selected sample quantiles in each cluster of observations.
#' 
#' @section Step 2. Three Types of Index Biomarker (??):
#' 
#' ## Optimal Dichotomizing Predictor(s) Selection via Dichotomizing Split Sample [optimSplit_dichotom]
#' 
#' Function [optimSplit_dichotom] identifies 
#' the optimal dichotomizing predictors using repeated sample splits on the *training set*.
#' 
#' Function [predict.optimSplit_dichotom] .. testing set
#' 
#' ## Functional Regression Indices [FRidx]
#' 
#' Function [FRidx] with option `nonlinear = FALSE` (default) .. training set
#' 
#' Function [predict.FRidx] .. testing set
#' 
#' Function [FRidx] with option `nonlinear = TRUE` .. training set
#' 
#' Function [predict.FRidx] .. testing set
#' 
#' @section Step 3. Bootstrap-Based Optimism Correction [BBC_dichotom]:
#' 
#' Function [BBC_dichotom]: Bootstrap-based optimism correction for dichotomizing selected predictor(s).
#' 
#' @example inst/extexample/exa_step1.R
#' @example inst/extexample/exa_step2.R
#' @example inst/extexample/exa_step3.R
#'
#' @references 
#' Selection of optimal quantile protein biomarkers based on cell-level immunohistochemistry data.
#' Misung Yi, Tingting Zhan, Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva. BMC Bioinformatics, 2023. \doi{10.1186/s12859-023-05408-8}
#' 
#' Quantile index biomarkers based on single-cell expression data.
#' Misung Yi, Tingting Zhan, Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva. 
#' Laboratory Investigation, 2023. \doi{10.1016/j.labinv.2023.100158}
#' 
'_PACKAGE'



