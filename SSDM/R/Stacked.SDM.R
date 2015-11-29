#' @include Ensemble.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent SSDMs
#'
#'This is an S4 class to represent SSDMs that assembles multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adaptive splines, generalized boosted regression model, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines) built for multiple species. It is obtained with
#'\code{\link{stack_modelling}} or \code{\link{stacking}}.
#'
#'@slot name character. Name of the SSDM (by default 'Species.SSDM').
#'@slot diversity.map raster. Local species richness map produced by the SSDM.
#'@slot uncertainty raster. Between-algorithm variance map.
#'@slot evaluation data frame. Evaluation of the SSDM (AUC, Kappa, omission
#'  rate, sensitivity, specificity, proportion of correctly predicted
#'  occurrences).
#'@slot variable.importance data frame. Relative importance of each variable in
#'  the SSDM.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'@slot enms list. List of ensemble SDMs used in the SSDM.
#'@slot parameters data frame. Parameters used to build the SSDM.
#'@slot algorithm.evaluation data frame. Evaluation of the algorithm averaging
#'  the metrics of all SDMs (AUC, Kappa, omission rate, sensitivity,
#'  specificity, proportion of correctly predicted occurrences).
#'
#'@seealso \linkS4class{Ensemble.SDM} an S4 class to represent ensemble SDMs,
#'  and \linkS4class{Algorithm.SDM} an S4 class to represent SDMs.
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
