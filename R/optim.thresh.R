# function SDMTools::optim.thresh() removed from CRAN as not maintained

.optim.thresh = function(obs,pred,threshold=101) {
  #input checks
  if (length(obs)!=length(pred)) stop('this requires the same number of observed & predicted values')
  if (length(threshold)==1 & threshold[1]>1) { thresholds = seq(0,1,length=threshold) } else { stop('inappropriate threshold value. See help file') }

  #deal with NAs
  if (length(which(is.na(c(obs,pred))))>0) {
    na = union(which(is.na(obs)),which(is.na(pred)))
    warning(length(na),' data points removed due to missing data')
    obs = obs[-na]; pred = pred[-na]
  }

  #define the n's and do checks
  n = length(obs); if (length(which(obs %in% c(0,1)))!=n) stop('observed values must be 0 or 1') #ensure observed are values 0 or 1

  #calculate the accuracy information for the data
  ac = .accuracy(obs,pred,threshold)

  #define the output list
  out = list()

  #calculate each of the optimum thresholds

  #min prediction of occurrences
  out$min.occurence.prediction = min(pred[which(obs==1)])

  #mean occurance prediction
  out$mean.occurence.prediction = mean(pred[which(obs==1)])

  #10% omission rate
  t.thresh = ac$threshold[which(rank(abs(0.1 - ac$omission),ties.method="min")==1)] #get the differences from calculated omission rates to 10%
  if (length(t.thresh)>1) { out$'10.percent.omission' = range(t.thresh) } else { out$'10.percent.omission' = t.thresh }

  #sensitivity = specificity
  t.thresh = ac$threshold[which(rank(abs(ac$sensitivity - ac$specificity),ties.method="min")==1)]
  if (length(t.thresh)>1) { out$'sensitivity=specificity' = range(t.thresh) } else { out$'sensitivity=specificity' = t.thresh }

  #max sensitivity + specificity
  t.thresh = ac$threshold[which(rank(-(ac$sensitivity + ac$specificity),ties.method="min")==1)]
  if (length(t.thresh)>1) { out$'max.sensitivity+specificity' = range(t.thresh) } else { out$'max.sensitivity+specificity' = t.thresh }

  #max Kappa
  t.thresh = ac$threshold[which(rank(-(ac$Kappa),ties.method="min")==1)]
  if (length(t.thresh)>1) { out$maxKappa = range(t.thresh) } else { out$maxKappa = t.thresh }

  #Max prop correct
  t.thresh = ac$threshold[which(rank(-(ac$prop.correct),ties.method="min")==1)]
  if (length(t.thresh)>1) { out$max.prop.correct = range(t.thresh) } else { out$max.prop.correct = t.thresh }

  #min ROC plot distance
  t.thresh = ac$threshold[which(rank((1-ac$sensitivity)^2 + (ac$specificity-1)^2,ties.method="min")==1)]
  if (length(t.thresh)>1) { out$min.ROC.plot.distance = range(t.thresh) } else { out$min.ROC.plot.distance = t.thresh }

  #return the values
  return(out)
}
