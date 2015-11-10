#' @include Ensemble.Niche.Model.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent a stack species distribution ensemble model of
#'multiple algorithms and multiple species
#'
#'This is an S4 class to represent a stack species distribution ensemble model
#'of multiple algorithms (among generalized linear model, general additive
#'model, multivariate adpatative splines, generalized boosted regression models,
#'classification tree analysis, random forest, maximum entropy, artificial
#'neural network, and support vector machines). It can be obtain with
#'\code{\link{Stack.Modelling}} or \code{\link{stacking}}.
#'
#'@slot name character. Name of the model (by default
#'  Stacked.Species.Niche.Model)
#'@slot diversity.map raster. Local species richness map of the model
#'@slot uncertainity raster. Intermodel variance map
#'@slot evaluation data frame. Evaluation of the model (threshold, AUC, omission
#'  rate, sensitivity, specificity, correct proportion and Kappa)
#'@slot variables.importance data frame. Relative percentage of importance for
#'  each variable used in the model
#'@slot algorithm.correlation data.frame. Correlation matrix inter algorithms
#'  calculated with pearson coefficient on habitat suitability map
#'@slot enms list. List with each ensemble model used for the stack model
#'@slot parameters data frame. Parameters used to realized the model
#' @slot algorithm.evaluation data frame. Evaluation of each algorihtm with the
#'  same metrics as evaluation slot
#'
#'@seealso \linkS4class{Ensemble.Niche.Model} an S4 class for ensemble models,
#'  and \linkS4class{Algorithm.Niche.Model} an S4 class for algorithm models.
#'
#' @export
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

# Class Generator
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
