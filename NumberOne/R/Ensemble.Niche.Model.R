#' @include Algorithm.Niche.Model.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent a specie distribution ensemble model of multiple
#'algorithms
#'
#'This is an S4 class to represent a specie distribution ensemble model of
#'multiple algorithms (among generalized linear model, general additive model,
#'multivariate adpatative splines, generalized boosted regression models,
#'classification tree analysis, random forest, maximum entropy, artificial
#'neural network, and support vector machines). It can be obtain with
#'\code{\link{Ensemble.Modelling}} or \code{\link{ensemble}}.
#'
#' @slot uncertainity raster. Intermodel variance map
#' @slot algorithm.correlation data.frame. Correlation matrix inter algorithms
#' @slot algorithm.evaluation
#'  calculated with pearson coefficient on habitat suitability map
#'
#' @seealso \linkS4class{Algorithm.Niche.Model} an S4 class for specie models
#'  with only one algorithm, and
#'  \linkS4class{Stack.Species.Ensemble.Niche.Model} an S4 class for stack
#'  species enemble models.
#'
#' @export
setClass('Ensemble.Niche.Model',
         contains = 'Niche.Model',
         representation(uncertainity = 'Raster',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame'),
         prototype(uncertainity = raster(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame()))

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
