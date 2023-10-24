test_that('stack modelling function', {
  data(Env)
  data(Occurrences)
  # parrallel
  i <- "species" # 'algorithms' and 'replicates' missing
  SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
                         Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                         Spcol = 'SPECIES', ensemble.thresh = 0,
                         verbose = FALSE,
                         cores = 2, parmode = i , tmp = FALSE)
  expect_is(SSDM, 'Stacked.SDM')
  # projection
  Env_new <- stack(Env[[1]]-0.3,Env[[2]],Env[[3]])
  SSDM_proj <- project(SSDM,Env_new,output.format='model')
  expect_is(SSDM_proj, 'Stacked.SDM')
  expect_false(all(is.na(values(SSDM_proj@diversity.map))))
})
