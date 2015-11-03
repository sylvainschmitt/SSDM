#' @include Ensemble.Niche.Model.R Stack.Niche.Model.R
#' @importFrom raster writeRaster
NULL

#' Save models
#'
#' Allow to save S4 \linkS4class{Ensemble.Niche.Model} and
#' \linkS4class{Stack.Species.Ensemble.Niche.Model} classes objects.
#'
#' @param enm Ensemble.Niche.Model. Ensemble model to be saved.
#' @param stack Stack.Species.Ensemble.Niche.Model.Stack species ensemble models
#'   to be saved.
#' @param name character. Folder name of the saved model.
#' @param directory character. Path to the directory containing the saved model
#'   folder, by default the current directory.
#'
#' @return Nothing in R environment. Create folders and tables and rasters
#'   associated to the model. Tables are in .csv and rasters in .grd/.gri.
#'
#' @seealso \code{\link{load.model}}
#'
#' @name save.model
#'
NULL

#' @rdname save.model
#' @export
setMethod('save.enm', 'Ensemble.Niche.Model', function (enm,
                                                        name = strsplit(enm@name, '.', fixed = T)[[1]][1],
                                                        directory = getwd()) {

  cat('Saving ensemble model results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/", name))
  dir.create(path = paste0(directory, "/", name,"/Rasters"))
  dir.create(path = paste0(directory, "/", name, "/Tables"))

  # Raster saving
  cat('   rasters ...')
  writeRaster(enm@projection[[1]], paste0(directory, "/", name, '/Rasters/Probability'), 'GTiff', overwrite = T)
  writeRaster(enm@uncertainity, paste0(directory, "/", name, '/Rasters/Uncertainity'), 'GTiff', overwrite = T)
  cat('saved \n')

  # Tables saving
  cat('   tables ...')
  write.csv(enm@evaluation, paste0(directory, "/", name, '/Tables/ENMeval'))
  write.csv(enm@algorithm.evaluation, paste0(directory, "/", name, '/Tables/AlgoEval'))
  write.csv(enm@algorithm.correlation, paste0(directory, "/", name, '/Tables/AlgoCorr'))
  write.csv(enm@variables.importance, paste0(directory, "/", name, '/Tables/VarImp'))
  write.csv(enm@data, paste0(directory, "/", name, '/Tables/Data'))
  write.csv(enm@name, paste0(directory, "/", name, '/Tables/Name'))
  write.csv(enm@parameters, paste0(directory, "/", name, '/Tables/Parameters'))
  cat('saved \n \n')
})

#' @rdname save.model
#' @export
setMethod('save.stack', 'Stack.Species.Ensemble.Niche.Model', function (stack, name = 'Stack', directory = getwd()) {

  cat('Saving stack species model results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/", name))
  directory = paste0(directory, "/", name)
  dir.create(path = paste0(directory, "/", 'Results'))
  dir.create(path = paste0(directory, "/", 'Results',"/Rasters"))
  dir.create(path = paste0(directory, "/", 'Results', "/Tables"))

  # Raster saving
  cat('   rasters ...')
  writeRaster(stack@diversity.map, paste0(directory, "/", 'Results', '/Rasters/Diversity'), 'GTiff', overwrite = T)
  writeRaster(stack@uncertainity, paste0(directory, "/", 'Results', '/Rasters/Uncertainity'), 'GTiff', overwrite = T)
  cat('saved \n')

  # Tables saving
  cat('   tables ...')
  write.csv(stack@evaluation, paste0(directory, "/", 'Results', '/Tables/StackEval'))
  write.csv(stack@algorithm.evaluation, paste0(directory, "/", 'Results', '/Tables/AlgoEval'))
  write.csv(stack@algorithm.correlation, paste0(directory, "/", 'Results', '/Tables/AlgoCorr'))
  write.csv(stack@variables.importance, paste0(directory, "/", 'Results', '/Tables/VarImp'))
  write.csv(stack@name, paste0(directory, "/", 'Results', '/Tables/Name'))
  write.csv(stack@parameters, paste0(directory, "/", 'Results', '/Tables/Parameters'))
  cat('saved \n\n')

  # ENMS saving
  cat('   enms ... \n\n')
  for (i in 1:length(stack@enms)) {save.enm(stack@enms[[i]], directory = directory)}
  cat('saved \n \n')
})
