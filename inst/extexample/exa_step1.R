### Data Preparation

library(survival)
data(Ki67, package = 'Qindex.data')
Ki67c = within(Ki67[complete.cases(Ki67), , drop = FALSE], expr = {
  marker = log1p(Marker); Marker = NULL
  PFS = Surv(RECFREESURV_MO, RECURRENCE)
})
(npt = length(unique(Ki67c$PATIENT_ID))) # 592

### Step 1: Cluster-Specific Sample Quantiles

Ki67q = clusterQp(marker ~ . - tissueID - inner_x - inner_y | PATIENT_ID, data = Ki67c)
stopifnot(is.matrix(Ki67q$marker))
head(Ki67q$marker, n = c(4L, 6L))

set.seed(234); id = sort.int(sample.int(n = npt, size = 480L))
Ki67q_0 = Ki67q[id, , drop = FALSE] # training set
Ki67q_1 = Ki67q[-id, , drop = FALSE] # test set
