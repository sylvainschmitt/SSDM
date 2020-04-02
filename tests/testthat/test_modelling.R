test_that('modelling function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  available.algo <- available.algo <- c("GLM", "GAM", "MARS", "GBM", "CTA", "RF", "ANN", "SVM")
  for(i in available.algo){
    cat(paste("testing", i, "...\n"))
    show_failure(SDM <- modelling(i, Occurrences, Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = F))
    expect_is(SDM, paste0(i,'.SDM'))
  }
})
