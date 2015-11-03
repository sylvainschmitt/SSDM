#' @include Algorithm.Niche.Model.R
#' @importFrom raster raster stack
NULL

setClass('Ensemble.Niche.Model',
         contains = 'Niche.Model',
         representation(uncertainity = 'Raster',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame'),
         prototype(uncertainity = raster(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame()))

#' @export
Ensemble.Niche.Model <- function(name = character(),
                                 projection = raster(),
                                 evaluation = data.frame(),
                                 variables.importance = data.frame(),
                                 data = data.frame(),
                                 uncertainity = raster(),
                                 algorithm.correlation = data.frame(),
                                 algorithm.evaluation = data.frame(),
                                 parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Ensemble.Niche.Model',
             name = name,
             projection = projection,
             evaluation = evaluation,
             variables.importance = variables.importance,
             data = data,
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             parameters = parameters))}
