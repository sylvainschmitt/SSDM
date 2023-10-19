test_that('project Ensemble.SDM', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  ESDM <- ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 2, Xcol = 'LONGITUDE', Ycol = 'LATITUDE', ensemble.thresh = 0, verbose = FALSE)
  Env_new <- stack(Env[[1]]-1,Env[[2]],Env[[3]])

  # cat(paste("test sequential ESDM projection \n"))
  ESDM_proj <- project(ESDM,Env_new,output.format='model')
  expect_is(ESDM_proj, 'Ensemble.SDM')
  expect_false(all(is.na(values(ESDM_proj@projection))))

  # cat(paste("test parallel ESDM projection \n"))
  ESDM_proj <- project(ESDM,Env_new,output.format='model',cores=2,minimal.memory=TRUE)
  expect_is(ESDM_proj, 'Ensemble.SDM')
  expect_false(all(is.na(values(ESDM_proj@projection))))

  # cat(paste("test ESDM projection with raster output \n"))
  ESDM_proj <- project(ESDM,Env_new,output.format='rasters')
  expect_is(ESDM_proj, 'list')
  expect_is(ESDM_proj$projection, 'RasterLayer')
  expect_false(all(is.na(values(ESDM_proj$projection))))
})
