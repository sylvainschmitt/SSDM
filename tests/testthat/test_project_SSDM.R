test_that('project Stacked.SDM', {
  data(Env)
  data(Occurrences)
  SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 2,
                          Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                          Spcol = 'SPECIES', ensemble.thresh = 0,
                          verbose = FALSE, cores = 0)
  
  Env_new <- stack(Env[[1]]-0.3,Env[[2]],Env[[3]])
  cat(paste("test sequential SSDM projection \n"))
  SSDM_proj <- project(SSDM,Env_new,output.format='model')
  expect_is(SSDM_proj, 'Stacked.SDM')
  expect_false(all(is.na(values(SSDM_proj@diversity.map))))
  
  cat(paste("test parallel SSDM projection \n"))
  SSDM_proj <- project(SSDM,Env_new,output.format='model',cores=2, minimal.memory=TRUE)
  expect_is(SSDM_proj, 'Stacked.SDM')
  expect_false(all(is.na(values(SSDM_proj@diversity.map))))
  
  cat(paste("test SSDM projection with raster output \n"))
  SSDM_proj <- project(SSDM,Env_new,output.format='rasters')
  expect_is(SSDM_proj, 'list')
  expect_is(SSDM_proj$diversity.map, 'RasterLayer')
  expect_false(all(is.na(values(SSDM_proj$diversity.map))))
})
