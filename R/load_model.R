#' @include Ensemble.SDM.R Stacked.SDM.R
#' @importFrom shiny incProgress
#' @importFrom raster raster stack
NULL

#' Function to load ensemble SDMs and SSDMs
#'
#' Allows to load S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' objects saved with their respective save function.
#'
#' @param name character. Name of the folder that contains the model to be
#'   loaded.
#' @param path character. Path to the directory containing the model to be
#'   loaded, by default the path to the current directory.
#' @param GUI logical. Don't take that argument into account (parameter for the
#'   user interface).
#'
#' @return The corresponding SDM object.
#'
#' @seealso \code{\link{save.model}}
#'
#' @name load.model
#'
NULL

#' @rdname load.model
#' @export
load_enm = function (name, path = getwd()) {
  path = paste0(path, '/', name)
  a = try(read.csv(paste0(path,'/Tables/AlgoCorr.csv'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty ! \n')
    a = data.frame()
  }
  enm = Ensemble.SDM(name = as.character(read.csv(paste0(path,'/Tables/Name.csv'))[1,2]),
                             projection = raster(paste0(path,'/Rasters/Probability.tif')),
                             uncertainty = try(raster(paste0(path,'/Rasters/uncertainty.tif'))),
                             evaluation = read.csv(paste0(path,'/Tables/ENMeval.csv'), row.names = 1),
                             algorithm.evaluation  = read.csv(paste0(path,'/Tables/AlgoEval.csv'), row.names = 1),
                             algorithm.correlation = a,
                             data = read.csv(paste0(path,'/Tables/Data.csv'), row.names = 1),
                             variable.importance = read.csv(paste0(path,'/Tables/VarImp.csv'), row.names = 1),
                             parameters = read.csv(paste0(path,'/Tables/Parameters.csv'), row.names = 1, colClasses = "character"))
  return(enm)
}

#' @rdname load.model
#' @export
load_stack = function (name = 'Stack', path = getwd(), GUI = F) {
  path = paste0(path, '/', name)
  a = try(read.csv(paste0(path,'/Stack/Tables/AlgoCorr.csv'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty ! \n')
    a = data.frame()
  }
  stack = Stacked.SDM(name = as.character(read.csv(paste0(path,'/Stack/Tables/Name.csv'))[1,2]),
                      diversity.map = raster(paste0(path,'/Stack/Rasters/Diversity.tif')),
                      endemism.map = raster(paste0(path,'/Stack/Rasters/Endemism.tif')),
                      uncertainty = raster(paste0(path,'/Stack/Rasters/uncertainty.tif')),
                      evaluation = read.csv(paste0(path,'/Stack/Tables/StackEval.csv'), row.names = 1),
                      variable.importance = read.csv(paste0(path,'/Stack/Tables/VarImp.csv'), row.names = 1),
                      algorithm.correlation = a,
                      algorithm.evaluation = read.csv(paste0(path,'/Stack/Tables/AlgoEval.csv'), row.names = 1),
                      enms = list(),
                      parameters = read.csv(paste0(path,'/Stack/Tables/Parameters.csv'), row.names = 1, colClasses = "character"))
  enms = list.dirs(paste0(path, '/Species'), recursive = F, full.names = F)
  if(GUI) {incProgress(1/(length(enms)+1), detail = 'stack main results')}
  for (i in 1:length(enms)) {
    cat(enms[i])
    enm = load_enm(enms[i], path = paste0(path, '/Species'))
    stack@enms[[enm@name]] = enm
    if(GUI) {incProgress(1/(length(enms)+1), detail = enm@name)}
    }
  return(stack)
}
