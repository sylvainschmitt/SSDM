#' @include Ensemble.SDM.R Stacked.SDM.R
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
  a = try(read.csv(paste0(path,'/Tables/AlgoCorr'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty !')
    a = data.frame()
  }
  enm = Ensemble.SDM(name = as.character(read.csv(paste0(path,'/Tables/Name'))[1,2]),
                             projection = raster(paste0(path,'/Rasters/Probability.tif')),
                             uncertainty = try(raster(paste0(path,'/Rasters/uncertainty.tif'))),
                             evaluation = read.csv(paste0(path,'/Tables/ENMeval'), row.names = 1),
                             algorithm.evaluation  = read.csv(paste0(path,'/Tables/AlgoEval'), row.names = 1),
                             algorithm.correlation = a,
                             data = read.csv(paste0(path,'/Tables/Data'), row.names = 1),
                             variable.importance = read.csv(paste0(path,'/Tables/VarImp'), row.names = 1),
                             parameters = read.csv(paste0(path,'/Tables/Parameters'), row.names = 1, colClasses = "character"))
  return(enm)
}

#' @rdname load.model
#' @export
load_stack = function (name = 'Stack', path = getwd(), GUI = F) {
  path = paste0(path, '/', name)
  stack = Stacked.SDM(name = as.character(read.csv(paste0(path,'/Results/Tables/Name'))[1,2]),
                      diversity.map = raster(paste0(path,'/Results/Rasters/Diversity.tif')),
                      endemism.map = raster(paste0(path,'/Results/Rasters/Endemism.tif')),
                      uncertainty = raster(paste0(path,'/Results/Rasters/uncertainty.tif')),
                      evaluation = read.csv(paste0(path,'/Results/Tables/StackEval'), row.names = 1),
                      variable.importance = read.csv(paste0(path,'/Results/Tables/VarImp'), row.names = 1),
                      algorithm.correlation = read.csv(paste0(path,'/Results/Tables/AlgoCorr'), row.names = 1),
                      algorithm.evaluation = read.csv(paste0(path,'/Results/Tables/AlgoEval'), row.names = 1),
                      enms = list(),
                      parameters = read.csv(paste0(path,'/Results/Tables/Parameters'), row.names = 1, colClasses = "character"))
  enms = list.dirs(path, recursive = F, full.names = F)
  enms = enms[-which(enms == 'Results')]
  if(GUI) {incProgress(1/(length(enms)+1), detail = 'stack main results')}
  for (i in 1:length(enms)) {
    enm = load_enm(enms[i], path = path)
    stack@enms[[enm@name]] = enm
    if(GUI) {incProgress(1/(length(enms)+1), detail = enm@name)}
    }
  return(stack)
}
