
### Step 3 (after Step 1 & 2)

Ki67q_0a = within.data.frame(Ki67q_0, expr = {
  FR = std_IQR(fr) 
  nlFR = std_IQR(nlfr)
  optS = std_IQR(marker[,'0.27'])
})
Ki67q_1a = within.data.frame(Ki67q_1, expr = {
  FR = std_IQR(fr_1)
  nlFR = std_IQR(nlfr_1)
  optS = std_IQR(marker[,'0.27']) 
})
# `optS`: use the best quantile but discard the cutoff identified by [optimSplit_dichotom]
# all models below can also be used on training data `Ki67q_0a`
# naive use
summary(coxph(PFS ~ NodeSt + Tstage + FR, data = Ki67q_1a))
summary(coxph(PFS ~ NodeSt + Tstage + nlFR, data = Ki67q_1a))
summary(coxph(PFS ~ NodeSt + Tstage + optS, data = Ki67q_1a))
# set.seed if necessary
summary(BBC_dichotom(PFS ~ NodeSt + Tstage ~ FR, data = Ki67q_1a))
# `NodeSt`, `Tstage`: predctors to be used as-is
# `FR` to be dichotomized
# set.seed if necessary
summary(BBC_dichotom(PFS ~ NodeSt + Tstage ~ nlFR, data = Ki67q_1a))
# set.seed if necessary
summary(BBC_dichotom(PFS ~ NodeSt + Tstage ~ optS, data = Ki67q_1a))
