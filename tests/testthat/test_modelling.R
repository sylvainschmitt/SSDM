test_that('modelling function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  available.algo <- available.algo <- c("GLM", "GAM", "MARS", "GBM", "CTA", "RF", "ANN", "SVM")
  for(i in available.algo){
    show_failure(SDM <- modelling(i, Occurrences, Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = F))
    expect_is(SDM, paste0(i,'.SDM'))
  }
  # projection
  Env_new <- stack(Env[[1]]-1,Env[[2]],Env[[3]])
  SDM_proj <- project(SDM, Env_new, output.format='model')
  expect_is(SDM_proj, 'SVM.SDM')
  expect_false(all(is.na(values(SDM_proj@projection))))
})
