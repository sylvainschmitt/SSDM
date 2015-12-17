#' @include Algorithm.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent an ensemble SDM
#'
#'This is an S4 class to represent an ensemble SDM from multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adaptive splines, generalized boosted regression model, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines). This S4 class is obtained with
#'\code{\link{ensemble_modelling}} or \code{\link{ensemble}}.
#'
#'@slot uncertainty raster. Between-algorithm variance map.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'@slot algorithm.evaluation data frame. Evaluation of the ensemble SDM (available
#'  metrics include AUC, Kappa, sensitivity, specificity and proportion of
#'  correctly predicted occurrences) and identification of the optimal threshold
#'  to convert the habitat suitability map into a binary presence/absence map.
#'
#'@seealso \linkS4class{Algorithm.SDM} an S4 class to represent an SDM based on
#'  a single algorithm, and \linkS4class{Stacked.SDM} an S4 class for SSDMs.
#'
#'@export
setClass('Ensemble.SDM',
         contains = 'SDM',
         representation(uncertainty = 'Raster',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame'),
         prototype(uncertainty = raster(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame()))

Ensemble.SDM <- function(name = character(),
                                 projection = raster(),
                                 evaluation = data.frame(),
                                 variable.importance = data.frame(),
                                 data = data.frame(),
                                 uncertainty = raster(),
                                 algorithm.correlation = data.frame(),
                                 algorithm.evaluation = data.frame(),
                                 parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Ensemble.SDM',
             name = name,
             projection = projection,
             evaluation = evaluation,
             variable.importance = variable.importance,
             data = data,
             uncertainty = uncertainty,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             parameters = parameters))}
