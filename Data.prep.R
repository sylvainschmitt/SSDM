##### Libraries ##### ----
library(spThin)

##### Occurences treatment ##### ----
treat.occ = function(Occurences, Env, Xcol, Ycol, Spcol = NULL, GeoRes = T, reso = max(res(Env@layers[[1]]))) {
  if (is.null(Spcol)) {
    Occurences$Sp = 1
    Spcol = 'NULL'
  }
  # Geographical resampling
  if (GeoRes) {
    cat('Geographical resampling \n\n')
    thin.result = thin(Occurences, long.col = Xcol, lat.col = Ycol, spec.col = 'Sp', 
                      thin.par = reso, reps = 1, locs.thinned.list.return = T, 
                      write.files = F, write.log.file = F, verbose = F)
    deleted = {}
    occ.indices = c(1:length(row.names(Occurences)))
    res.indices = as.numeric(row.names(thin.result[[1]]))
    for (i in 1:length(occ.indices)) {if(!(occ.indices[i] %in% res.indices)) {deleted = c(deleted, occ.indices[i])}}
    Occurences = Occurences[-deleted,]
  }
  if (Spcol == 'NULL') {Occurences = Occurences[-which(names(Occurences)=='Sp')]}
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
    reso = res(Env@layers[[1]])
    # Coarser resolution measurement
    for (i in 1:length(Env@layers)){
      reso[1] = max((res(Env@layers[[i]]))[1],reso[1])
      reso[2] = max((res(Env@layers[[i]]))[1],reso[2])
    }
    # Refine all stack resolution
    for (i in 1:length(Env@layers)) {
      Env[[i]] = aggregate(Env[[i]], fact = c(res(Env[[i]])/reso[1],res(Env[[i]])/reso[2]), fun = max)
    }
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
