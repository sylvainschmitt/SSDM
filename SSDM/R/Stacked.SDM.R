#' @include Ensemble.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent SSDMs
#'
#'This is an S4 class to represent SSDMs of multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adpatative splines, generalized boosted regression models, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines). It can be obtain with \code{\link{stack_modelling}}
#'or \code{\link{stacking}}.
#'
#'@slot name character. Name of the SDM (by default Species.SSDMs).
#'@slot diversity.map raster. Local species richness map of the Stacked SDMs
#'@slot uncertainty raster. Between-algorithm variance map.
#'@slot evaluation data frame. Evaluation of the Stacked SDMs (threshold, AUC, omission
#'  rate, sensitivity, specificity, correct proportion and Kappa).
#'@slot variable.importance data frame. Relative percentage of importance for
#'  each variable used in the model.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'@slot enms list. List with each ensemble SDM used for the stacked SDMs.
#'@slot parameters data frame. Parameters used to realized the stacked SDMs.
#'@slot algorithm.evaluation data frame. Calculated with Pearson's coefficient
#'  from habitat suitability map.
#'
#'@seealso \linkS4class{Ensemble.SDM} an S4 class for ensemble SDMs, and
#'  \linkS4class{Algorithm.SDM} an S4 class for SDMs.
#'
#'@export
setClass('Stacked.SDM',
         representation(name = 'character',
                        diversity.map = 'Raster',
                        uncertainty = 'Raster',
                        evaluation = 'data.frame',
                        variable.importance = 'data.frame',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame',
                        enms = 'list',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   diversity.map = raster(),
                   uncertainty = raster(),
                   evaluation = data.frame(),
                   variable.importance = data.frame(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame(),
                   enms = list(),
                   parameters = data.frame()))

# Class Generator
Stacked.SDM <- function(name = character(),
                                               diversity.map = raster(),
                                               uncertainty = raster(),
                                               evaluation = data.frame(),
                                               variable.importance = data.frame(),
                                               algorithm.correlation = data.frame(),
                                               algorithm.evaluation = data.frame(),
                                               enms = list(),
                                               parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Stacked.SDM',
             name = name,
             diversity.map = diversity.map,
             evaluation = evaluation,
             variable.importance = variable.importance,
             uncertainty = uncertainty,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             enms = enms,
             parameters = parameters))}
