#' @include Ensemble.SDM.R Stacked.SDM.R
#' @importFrom raster raster stack
NULL

#' Function to load ensemble SDMs and SSDMs
#'
#' Allow to load S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' objects saved with their respective save function.
#'
#' @param name character. Name of the folder that contains the model to be loaded.
#' @param directory character. Path to the directory containing the saved model
#'   folder, by default the current directory.
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
load_enm = function (name, directory = getwd()) {
  directory = paste0(directory, '/', name)
  a = try(read.csv(paste0(directory,'/Tables/AlgoCorr'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty !')
    a = data.frame()
  }
  enm = Ensemble.SDM(name = as.character(read.csv(paste0(directory,'/Tables/Name'))[1,2]),
                             projection = raster(paste0(directory,'/Rasters/Probability.tif')),
                             uncertainty = try(raster(paste0(directory,'/Rasters/uncertainty.tif'))),
                             evaluation = read.csv(paste0(directory,'/Tables/ENMeval'), row.names = 1),
                             algorithm.evaluation  = read.csv(paste0(directory,'/Tables/AlgoEval'), row.names = 1),
                             algorithm.correlation = a,
                             data = read.csv(paste0(directory,'/Tables/Data'), row.names = 1),
                             variable.importance = read.csv(paste0(directory,'/Tables/VarImp'), row.names = 1),
                             parameters = read.csv(paste0(directory,'/Tables/Parameters'), row.names = 1, colClasses = "character"))
  return(enm)
}

#' @rdname load.model
#' @export
load_stack = function (name = 'Stack', directory = getwd(), GUI = F) {
  directory = paste0(directory, '/', name)
  stack = Stacked.SDM(name = as.character(read.csv(paste0(directory,'/Results/Tables/Name'))[1,2]),
                                             diversity.map = raster(paste0(directory,'/Results/Rasters/Diversity.tif')),
                                             uncertainty = raster(paste0(directory,'/Results/Rasters/uncertainty.tif')),
                                             evaluation = read.csv(paste0(directory,'/Results/Tables/StackEval'), row.names = 1),
                                             variable.importance = read.csv(paste0(directory,'/Results/Tables/VarImp'), row.names = 1),
                                             algorithm.correlation = read.csv(paste0(directory,'/Results/Tables/AlgoCorr'), row.names = 1),
                                             algorithm.evaluation = read.csv(paste0(directory,'/Results/Tables/AlgoEval'), row.names = 1),
                                             enms = list(),
                                             parameters = read.csv(paste0(directory,'/Results/Tables/Parameters'), row.names = 1, colClasses = "character"))
  enms = list.dirs(directory, recursive = F, full.names = F)
  enms = enms[-which(enms == 'Results')]
  if(GUI) {incProgress(1/(length(enms)+1), detail = 'stack main results')}
  for (i in 1:length(enms)) {
    enm = load_enm(enms[i], directory = directory)
    stack@enms[[enm@name]] = enm
    if(GUI) {incProgress(1/(length(enms)+1), detail = enm@name)}
    }
  return(stack)
}
