##### Libraries ##### ----
library(spThin)

##### Occurences treatment ##### ----
load.occ = function(directory = getwd(), Env, file = NULL, ...,
                     Xcol = 'Longitude', Ycol = 'Latitude', Spcol = NULL, 
                     GeoRes = T, reso = max(res(Env@layers[[1]]))) {
  pdir = getwd()
  cat('Occurences loading \n')
  setwd(directory)
  if (is.null(file)) {file = as.character(list.files(pattern = '.csv')[[1]])}
  Occurences = read.csv2(file = file, ...)  # Occ = occurences
  
  # Geographical resampling
  if (is.null(Spcol)) {
    Occurences$Sp = 1
    Spcol = 'NULL'
  }
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
  
  setwd(pdir)
  return(Occurences)
}

##### Environment variables treatment ##### ----
# 2 -  Variables treatment #
# Without multiscaling for the moment
load.var <- function (directory = getwd(), files = NULL, 
                      format = c('.grd','.tif','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img'),
                      factors = NULL, Norm = T, tmp = T) {
  pdir = getwd()
  cat('Variables loading \n')
  setwd(directory)
  Env = stack()
  
  # Rasters loading
  files.null = files
  if (is.null(files)) {files.null = T} else {files.null = F}
  for (j in 1:length(format)) {
    if(files.null) {
      files = list.files(pattern = paste0('.',format[j],'$'))
    }
    if (length(files) > 0) {
      for (i in 1:length(files)){
        Raster = raster(files[[i]])
        # Extent and resolution check
        reso = res(Raster)
        extent = extent(Raster)
        if (j == 1  && i == 1) {
          resostack = reso
          extentstack = extent
        } else {
          resostack[1] = max(reso[1], resostack[1])
          resostack[2] = max(reso[2], resostack[2])
          # Extent and resolution adpatation
          extentstack@xmin = max(extentstack@xmin, extent@xmin)
          extentstack@xmax = min(extentstack@xmax, extent@xmax)
          extentstack@ymin = max(extentstack@ymin, extent@ymin)
          extentstack@ymax = min(extentstack@ymax, extent@ymax)
        }
      }
    }
  }
  
  cat('Variables treatment \n')
  
  cat('   resolution and extent adaptation...')
  for (j in 1:length(format)) {
    if(files.null) {
      files = list.files(pattern = paste0('.',format[j],'$'))
    }
    if (length(files) > 0) {
      for (i in 1:length(files)){
        Raster = raster(files[[i]])
        Raster = reclassify(Raster, c(-Inf,-900,NA))
        Raster = crop(Raster, extentstack)
        names(Raster) = as.character(strsplit(files[i],paste0('.',format[j]))[1])
        if (names(Raster) %in% factors) {
          Raster = as.factor(Raster)
          row.names(Raster@data@attributes[[1]]) = Raster@data@attributes[[1]]$ID
          fun = max
        } else {fun = mean}
        if (round(res(Raster)[1], digits = 6) != round(resostack[1], digits = 6) || round(res(Raster)[2], digits = 6) != round(resostack[2], digits = 6)) {
          cat(c((res(Raster)[1]/resostack[1]),(res(Raster)[2]/resostack[2])))
          Raster = aggregate(Raster, fact = (res(Raster)[1]/resostack[1]), fun = fun)
          }
        Env = stack(Env, Raster)
      }
    }
  }
  cat('   done... \n')
  
  
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
  
  setwd(pdir)
  
  # Temporary files
  if (tmp) {
    if (!("./.rasters" %in% list.dirs())) (dir.create('./.rasters'))
    setwd(".rasters")
    for (i in 1:length(Env@layers)) {Env[[i]] = writeRaster(Env[[i]], names(Env[[i]]), overwrite = T)}
    setwd(pdir)
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
