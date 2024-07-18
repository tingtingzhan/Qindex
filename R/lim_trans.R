

lim_trans <- function(old, new) {
  max_old <- max(old, na.rm = TRUE)
  min_old <- min(old, na.rm = TRUE)
  max_new <- max(new, na.rm = TRUE)
  min_new <- min(new, na.rm = TRUE)
  
  slope <- (max_old - min_old) / (max_new - min_new)
  intercept <- min_old - min_new * slope
  
  stopifnot(all.equal((range(old, na.rm = TRUE) - intercept) / slope, range(new, na.rm = TRUE)))
  
  return(c(intercept = intercept, slope = slope))
  
}