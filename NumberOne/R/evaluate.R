#' @include Algorithm.Niche.Model.R
#' @importFrom raster stack extract
#' @importFrom SDMTools optim.thresh accuracy
NULL

#'Method to evaluate a model
#'
#'This is a method to evaluate models from S4
#'\linkS4class{Algorithm.Niche.Model} or \linkS4class{Ensemble.Niche.Model}
#'classes objects.
#'
#'@param obj Algorithm.Niche.Model or Ensemble.Niche.Model. Model to evaluate
#'@param metric character. Method used to calculate the optimal threshold to
#'  compute the confusion matrix used for the evaluation of the model (see
#'  details below)
#'@param thresh numeric. Binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'
#'@details \strong{Metric :} choice of the method used to compute binary map
#'threshold and confusion matrix : \describe{ \item{Kappa}{maximizes the model
#'Kappa value} \item{TSS}{\strong{True Skill Statistic} maximizes the
#'sensitivity and specificity sum} \item{CCR}{maximizes the correct predicted
#'observations proportion} \item{SES}{using the sensitivty specificity
#'equality} \item{LW}{using the lowest occurence prediction probability}
#'\item{ROC}{minimizing the distance between the ROC plot (receiving operative
#'curve) and the upper left coin (1,1)} }
#'
#'@return S4 \linkS4class{Algorithm.Niche.Model} or
#'  \linkS4class{Ensemble.Niche.Model} classes objects. viewable with
#'  \code{\link{plot.model}} method depending on the input.
#'
#' @examples
#'\dontrun{
#' evaluate(GLM)
#'}
#'
#'@seealso \code{\link{Modelling}} for specie distribution modelling with one
#'  algorithm, \code{\link{Ensemble.Modelling}} for specie ensemble modelling
#'  with multiple algorithms
#'
#'@export
setMethod("evaluate", "Niche.Model", function(obj, thresh = 1001, metric = 'SES') {
  data = obj@data[which(!obj@data$Train),]
  predicted.values = extract(obj@projection, data[c('X','Y')])
  predicted.values[which(is.na(predicted.values))] = 0
  # Threshold computing
  metric = switch(metric,
                  'Kappa' = 'maxKappa',
                  'CCR' = 'max.prop.correct',
                  'TSS' = 'max.sensitivity+specificity',
                  'SES' = 'sensitivity=specificity',
                  'LW' = 'min.occurence.prediction',
                  'ROC' = 'min.ROC.plot.distance')
  threshold = optim.thresh(data$Presence, predicted.values, thresh)
  threshold = mean(threshold[[which(names(threshold) == metric)]])

  obj@evaluation = accuracy(data$Presence, predicted.values, threshold)
  row.names(obj@evaluation) = "Evaluation"
  return(obj)})
