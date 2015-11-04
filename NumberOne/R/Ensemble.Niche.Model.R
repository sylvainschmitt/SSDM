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
#'\code{\link{Ensemble.Modelling}} or \code{\link{ensemble,Algorithm.Niche.Model-method}}.
#'
#'@slot name character. Name of the model (by default
#'  Specie.Ensemble.Niche.Model)
#'@slot projection raster. Habitat suitability map of the model
#'@slot uncertainity raster. Intermodel variance map
#'@slot evaluation data frame. Evaluation of the model (threshold, AUC, omission
#'  rate, sensitivity, specificity, correct proportion and Kappa)
#'@slot algorithm.evluation data frame. Evaluation of each algorihtm with the
#'  same metrics as evaluation slot
#'@slot variables.importance data frame. Relative percentage of importance for
#'  each variable used in the model
#'@slot algorithm.correlation data.frame. Correlation matrix inter algorithms
#'  calculated with pearson coefficient on habitat suitability map
#'@slot data data frame. Data used to realized the model
#'@slot parameters data frame. Parameters used to realized the model
#'
#'@seealso \linkS4class{Algorithm.Niche.Model} an S4 class for specie models
#'  with only one algorithm, and
#'  \linkS4class{Stack.Species.Ensemble.Niche.Model} an S4 class for stack
#'  species enemble models.
#'
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
