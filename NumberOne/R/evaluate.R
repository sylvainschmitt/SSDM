#' @include Algorithm.Niche.Model.R
#' @importFrom raster stack extract
#' @importFrom SDMTools optim.thresh accuracy
NULL

#' @export
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
