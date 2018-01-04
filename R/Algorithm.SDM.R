#' @include SDM.R
#' @import methods
NULL

#' An S4 class to represent an SDM based on a single algorithm
#'
#' This is an S4 class to represent an SDM based on a single algorithm (including
#' generalized linear model, general additive model, multivariate adpative
#' splines, generalized boosted regression model, classification tree analysis,
#' random forest, maximum entropy, artificial neural network, and support vector
#' machines). This S4 class is obtained with \code{\link{modelling}}.
#'
#' @slot name character. Name of the SDM (by default Species.SDM).
#' @slot projection raster. Habitat suitability map produced by the SDM.
#' @slot binary raster. Presence/Absence binary map produced by the SDM.
#' @slot evaluation data frame. Evaluation of the SDM (available metrics include
#'  AUC, Kappa, sensitivity, specificity and proportion of correctly predicted
#'  occurrences) and identification of the optimal threshold to convert the
#'  habitat suitability map into a binary presence/absence map.
#' @slot variable.importance data frame. Relative importance of
#'  each variable in the SDM.
#' @slot data data frame. Data used to build the SDM.
#' @slot parameters data frame. Parameters used to build the SDM.
#'
#' @seealso \linkS4class{Ensemble.SDM} an S4 class for ensemble SDMs,
#'  and \linkS4class{Stacked.SDM} an S4 class for SSDMs.
#'
#' @export
setClass('Algorithm.SDM',
         contains = 'SDM')

# Class generator
Algorithm.SDM <- function(algorithm = 'Algorithm',
                          name = character(),
                          projection = raster(),
                          binary = raster(),
                          evaluation = data.frame(),
                          variable.importance = data.frame(),
                          data = data.frame(),
                          parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  object.class <- paste0(algorithm,'.SDM')
  return(new(object.class,
             name = name,
             binary = binary,
             projection = projection,
             evaluation = evaluation,
             variable.importance = variable.importance,
             data = data,
             parameters = parameters))
}

setClass('GLM.SDM',
         contains = 'Algorithm.SDM')

setClass('GAM.SDM',
         contains = 'Algorithm.SDM')

setClass('MARS.SDM',
         contains = 'Algorithm.SDM')

setClass('CTA.SDM',
         contains = 'Algorithm.SDM')

setClass('GBM.SDM',
         contains = 'Algorithm.SDM')

setClass('RF.SDM',
         contains = 'Algorithm.SDM')

setClass('MAXENT.SDM',
         contains = 'Algorithm.SDM')

setClass('ANN.SDM',
         contains = 'Algorithm.SDM')

setClass('SVM.SDM',
         contains = 'Algorithm.SDM')
