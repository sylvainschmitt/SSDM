##### Libraries ##### ----
library(spThin)

##### Occurences treatment ##### ----
treat.occ = function(Occurences, Env, Xcol, Ycol, Spcol, GeoRes = T, reso = max(res(Env@layers[[1]]))) {
  # Geographical resampling
  if (GeoRes) {
    cat('Geographical resampling \n\n')
    thin.result = thin(Occurences, long.col = Xcol, lat.col = Ycol, spec.col = Spcol, 
                      thin.par = reso, reps = 1, locs.thinned.list.return = T, 
                      write.files = F, write.log.file = F, verbose = F)
    deleted = {}
    occ.indices = c(1:length(row.names(Occurences)))
    res.indices = as.numeric(row.names(thin.result[[1]]))
    for (i in 1:length(occ.indices)) {if(!(occ.indices[i] %in% res.indices)) {deleted = c(deleted, occ.indices[i])}}
    Occurences = Occurences[-deleted,]
  }
  return(Occurences)
}

##### Environment variables treatment ##### ----
# 2 -  Variables treatment #
# Without multiscaling for the moment
treat.var <- function (Env, Reso.adapt = T, Norm = T) {
  cat('Variables treatment \n')
  
  # Resolution
  if(Reso.adapt) {
    cat('   resolution adaptation \n')
    reso = max(res(Env@layers[[1]]))
    # Coarser resolution measurement
    for (i in 1:length(Env@layers)){
      reso = max((res(Env@layers[[i]])),reso)
    }
    # Refine all stack resolution
    res(Env) = reso
  }
  
  # Normalizing variables
  if (Norm) {
    cat('   normalizing continuous variables \n\n')
    for (i in 1:length(Env@layers)) {
      #For not categorical variables
      if (!Env[[i]]@data@isfactor) { 
        Env[[i]] = Env[[i]]/Env[[i]]@data@max
      }
    }
  }
  
  return(Env)
}

##### Mask creation ##### ----
env.mask = function(Env) {
  mask = reclassify(Env[[1]], c(-Inf,Inf,0))
  for (i in 1:length(Env@layers)) {
    mask = mask + reclassify(Env[[i]], c(-Inf,Inf,1))
  }
  mask = mask / length(Env@layers)
  mask = reclassify(mask, c(-Inf,1,NA), right = F)
  return(mask)
}
