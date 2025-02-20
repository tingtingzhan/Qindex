

.onAttach <- function(libname, pkgname = 'Qindex') {

  # NOTHING NEEDED HERE!!!
  
  # use Depends
  # and 
  # force @import
  
  # using base::library gives me error on devtools::check
  
  # base::requireNamespace
  # .. requires end-user to library
  # .. prints message 'Loading required namespace: pkg'
  
  # base::require
  
  #require(spatstat.grouped)
  #require(maxEff) # previous [optim_split_etc].
  #require(boc) # previous [BBC_dichotom] function 
  #require(gam.matrix) # previous [Qindex] function
  
}
