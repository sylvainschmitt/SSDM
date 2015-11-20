#' @include Ensemble.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent SSDMs
#'
#'This is an S4 class to represent SSDMs that assembles multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adaptative splines, generalized boosted regression model, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines). It is obtained with \code{\link{stack_modelling}}
#'or \code{\link{stacking}}.
#'
#'@slot name character. Name of the SDM (by default 'Species.SSDM').
#'@slot diversity.map raster. Local species richness map of the SSDM.
#'@slot uncertainty raster. Between-algorithm variance map.
#'@slot evaluation data frame. Evaluation of the SSDM (AUC, omission rate,
#'  sensitivity, specificity, correct proportion and Kappa).
#'@slot variable.importance data frame. Relative importance of each
#'  variable in the SSDM.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'@slot enms list. List with each ensemble SDM used for the SSDM.
#'@slot parameters data frame. Parameters used to realized the SSDM.
#'@slot algorithm.evaluation data frame. Evaluation of the algorithm by is mean
#'  performance among all SDMs (AUC, omission rate, sensitivity, specificity,
#'  correct proportion and Kappa).
#'
#'@seealso \linkS4class{Ensemble.SDM} an S4 class to represent ensemble SDMs, and an
#'  \linkS4class{Algorithm.SDM} an S4 class to represent SDMs.
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
