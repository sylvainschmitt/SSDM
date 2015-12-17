#' @include Ensemble.SDM.R Stacked.SDM.R checkargs.R
#' @importFrom raster writeRaster
NULL

#' Save ensemble SDMs and SSDMs
#'
#' Allows to save S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' class objects.
#'
#' @param enm Ensemble.SDM. Ensemble SDM to be saved.
#' @param stack Stacked.SDM. SSDM to be saved.
#' @param name character. Folder name of the model to save.
#' @param path character. Path to the directory chosen to save the SDM,
#'   by default the path to the current directory.
#' @param verbose logical. If set to true, allows the function to print text in the
#'   console.
#' @param GUI logical. Don't take that argument into account (parameter for the
#'   user interface).
#'
#' @return Nothing in R environment. Creates folders, tables and rasters
#'   associated to the SDM. Tables are in .csv and rasters in .grd/.gri.
#'
#' @seealso \code{\link{load.model}}
#'
#' @name save.model
#'
NULL

#' @rdname save.model
#' @export
setGeneric('save.enm', function (enm, name = strsplit(enm@name, '.', fixed = T)[[1]][1],
                                 path = getwd(), verbose = T, GUI = F) {return(standardGeneric('save.enm'))})

#' @rdname save.model
#' @export
setMethod('save.enm', 'Ensemble.SDM', function (enm,
                                                        name = strsplit(enm@name, '.', fixed = T)[[1]][1],
                                                        path = getwd(),
                                                        verbose = T, GUI = F) {
  # Check arguments
  .checkargs(enm = enm, name = name, path = path,  verbose = verbose, GUI = GUI)

  if(verbose){cat('Saving ensemble model results \n')}
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  dir.create(path = paste0(path, "/", name,"/Rasters"))
  dir.create(path = paste0(path, "/", name, "/Tables"))

  # Raster saving
  if(verbose){cat('   rasters ...')}
  writeRaster(enm@projection[[1]], paste0(path, "/", name, '/Rasters/Probability'), 'GTiff', overwrite = T)
  writeRaster(enm@uncertainty, paste0(path, "/", name, '/Rasters/uncertainty'), 'GTiff', overwrite = T)
  if(verbose){cat('saved \n')}

  # Tables saving
  if(verbose){cat('   tables ...')}
  write.csv(enm@evaluation, paste0(path, "/", name, '/Tables/ENMeval'))
  write.csv(enm@algorithm.evaluation, paste0(path, "/", name, '/Tables/AlgoEval'))
  write.csv(enm@algorithm.correlation, paste0(path, "/", name, '/Tables/AlgoCorr'))
  write.csv(enm@variable.importance, paste0(path, "/", name, '/Tables/VarImp'))
  write.csv(enm@data, paste0(path, "/", name, '/Tables/Data'))
  write.csv(enm@name, paste0(path, "/", name, '/Tables/Name'))
  write.csv(enm@parameters, paste0(path, "/", name, '/Tables/Parameters'))
  if(verbose){cat('saved \n \n')}
})

#' @rdname save.model
#' @export
setGeneric('save.stack', function (stack, name = 'Stack', path = getwd(), verbose = T, GUI = F) {return(standardGeneric('save.stack'))})

#' @rdname save.model
#' @export
setMethod('save.stack', 'Stacked.SDM', function (stack, name = 'Stack',
                                                                        path = getwd(),
                                                                        verbose = T, GUI = F) {
  # Check arguments
  .checkargs(stack = stack, name = name, path = path,  verbose = verbose, GUI = GUI)

  if(verbose){cat('Saving stack species model results \n')}
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  path = paste0(path, "/", name)
  dir.create(path = paste0(path, "/", 'Results'))
  dir.create(path = paste0(path, "/", 'Results',"/Rasters"))
  dir.create(path = paste0(path, "/", 'Results', "/Tables"))

  # Raster saving
  if(verbose){cat('   rasters ...')}
  writeRaster(stack@diversity.map, paste0(path, "/", 'Results', '/Rasters/Diversity'), 'GTiff', overwrite = T)
  writeRaster(stack@endemism.map, paste0(path, "/", 'Results', '/Rasters/Endemism'), 'GTiff', overwrite = T)
  writeRaster(stack@uncertainty, paste0(path, "/", 'Results', '/Rasters/uncertainty'), 'GTiff', overwrite = T)
  cat('saved \n')

  # Tables saving
  if(verbose){cat('   tables ...')}
  write.csv(stack@evaluation, paste0(path, "/", 'Results', '/Tables/StackEval'))
  write.csv(stack@algorithm.evaluation, paste0(path, "/", 'Results', '/Tables/AlgoEval'))
  write.csv(stack@algorithm.correlation, paste0(path, "/", 'Results', '/Tables/AlgoCorr'))
  write.csv(stack@variable.importance, paste0(path, "/", 'Results', '/Tables/VarImp'))
  write.csv(stack@name, paste0(path, "/", 'Results', '/Tables/Name'))
  write.csv(stack@parameters, paste0(path, "/", 'Results', '/Tables/Parameters'))
  cat('saved \n\n')

  # ENMS saving
  if(verbose){cat('   enms ... \n\n')}
  for (i in 1:length(stack@enms)) {save.enm(stack@enms[[i]], path = path, verbose = verbose, GUI = GUI)}
  if(verbose){cat('saved \n \n')}
})