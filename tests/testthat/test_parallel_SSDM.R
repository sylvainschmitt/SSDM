test_that('Stack parallelization', {
  data(Env)
  data(Occurrences)
  for(i in c('species', 'algorithms', 'replicates')){
    SSDM <- stack_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 1,
                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE', Spcol = 'SPECIES',
                            ensemble.thresh = c(0), method="pSSDM", uncertainty = TRUE,
                            verbose = FALSE, cores = 2, parmode = i , minimal.memory = TRUE, tmp = FALSE)
    expect_is(SSDM, 'Stacked.SDM')
  }
})
