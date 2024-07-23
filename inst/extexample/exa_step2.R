
### Step 2 (after Step 1)

## Step 2a: Linear Functional Regression Indices
(fr = FRidx(PFS ~ marker, data = Ki67q_0))
stopifnot(all.equal.numeric(c(fr), predict(fr)))
integrandSurface(fr)
fr_1 = predict(fr, newdata = Ki67q_1)
integrandSurface(fr, newdata = Ki67q_1)

## Step 2b: Non-Linear Functional Regression Indices
(nlfr = FRidx(PFS ~ marker, data = Ki67q_0, nonlinear = TRUE))
stopifnot(all.equal.numeric(c(nlfr), predict(nlfr)))
integrandSurface(nlfr)
if (FALSE) {
  integrandSurface(nlfr, subj_vis = NULL)
  vis.gam(nlfr@gam)
  #autoplot.vis_gam_(vis_gam_(nlfr@gam))
} # these are mathematically identical
nlfr_1 = predict(nlfr, newdata = Ki67q_1)
integrandSurface(nlfr, newdata = Ki67q_1)

### Step 2c: Optimal Dichotomizing
set.seed(14837); (m1 = optimSplit_dichotom(
  PFS ~ marker, data = Ki67q_0, nsplit = 20L, 
  include = (highX > .15 & highX < .85), top = 2L)) 
predict(m1)
predict(m1, boolean = FALSE)
predict(m1, newdata = Ki67q_1) # using both quantile selected AND cutoff identified.
