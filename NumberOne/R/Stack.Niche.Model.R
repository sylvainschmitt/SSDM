#' @include Ensemble.Niche.Model.R
#' @importFrom raster raster stack
NULL

##### Stack Species Ensemble Niche Model Class ##### -----

# 1 - Class definition #
setClass('Stack.Species.Ensemble.Niche.Model',
         representation(name = 'character',
                        diversity.map = 'Raster',
                        uncertainity = 'Raster',
                        evaluation = 'data.frame',
                        variables.importance = 'data.frame',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame',
                        enms = 'list',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   diversity.map = raster(),
                   uncertainity = raster(),
                   evaluation = data.frame(),
                   variables.importance = data.frame(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame(),
                   enms = list(),
                   parameters = data.frame()))

# 2 - Class creation function #
#' @export
Stack.Species.Ensemble.Niche.Model <- function(name = character(),
                                               diversity.map = raster(),
                                               uncertainity = raster(),
                                               evaluation = data.frame(),
                                               variables.importance = data.frame(),
                                               algorithm.correlation = data.frame(),
                                               algorithm.evaluation = data.frame(),
                                               enms = list(),
                                               parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Stack.Species.Ensemble.Niche.Model',
             name = name,
             diversity.map = diversity.map,
             evaluation = evaluation,
             variables.importance = variables.importance,
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             enms = enms,
             parameters = parameters))}
