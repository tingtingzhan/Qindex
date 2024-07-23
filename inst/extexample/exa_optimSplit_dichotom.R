
library(survival)
data(pbc, package = 'survival') # see more details from ?survival::pbc
dim(pbc)
head(pbc2 <- within.data.frame(subset(pbc, status != 1L), expr = {
  death = (status == 2L)
  trt = structure(trt, levels = c('D-penicillmain', 'placebo'), class = 'factor')
  trt = relevel(trt, ref = 'placebo')
}))

sapply(pbc2, class)

nn <- nrow(pbc2)
set.seed(seed = 1151); id_train <- sample(seq_len(nn), size = floor(.8*nn))
dim(pbc2_train <- pbc2[id_train, ])
dim(pbc2_test <- pbc2[setdiff(seq_len(nn), id_train), ])

if (FALSE) {
  set.seed(seed = 13542); m0 = optimSplit_dichotom_VERYOLD(
    Surv(time, death) ~ bili + chol + albumin + copper + alk.phos + ast + trig + platelet + protime, 
    data = pbc2_train, nsplit = 20L, top = 2L) 
  head(m0, n = 10L)
  attr(m0, 'top')
}

set.seed(seed = 13542); (m1 = optimSplit_dichotom(
  Surv(time, death) ~ bili + chol + albumin + copper + alk.phos + ast + trig + platelet + protime, 
  data = pbc2_train, nsplit = 20L, top = 2L)) 
predict(m1)
predict(m1, data = pbc2_test)



