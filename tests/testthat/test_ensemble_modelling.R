test_that('ensemble modelling function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  # parrallel
  ESDM <- ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 1,
                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                            ensemble.thresh = c(0), verbose = FALSE,
                            cores = 2, tmp = TRUE)
  expect_is(ESDM, 'Ensemble.SDM')
  # projection
  Env_new <- stack(Env[[1]]-1,Env[[2]],Env[[3]])
  ESDM_proj <- project(ESDM,Env_new,output.format='model')
  expect_is(ESDM_proj, 'Ensemble.SDM')
  expect_false(all(is.na(values(ESDM_proj@projection))))
})
