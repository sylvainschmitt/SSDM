test_that('Ensemble parallelization', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  ESDM <- ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 10,
                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                            ensemble.thresh = c(0), verbose = FALSE, cores = 2,
                            minimal.memory = TRUE, tmp = TRUE)
  expect_is(ESDM, 'Ensemble.SDM')
})
