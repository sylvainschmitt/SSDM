#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom SDMTools optim.thresh accuracy
#' @importFrom stats aggregate.data.frame cor glm glm.control rbinom runif sd var
#' @importFrom utils lsf.str read.csv read.csv2 tail write.csv
NULL

setGeneric('evaluate', function(obj, ...) {return(standardGeneric('evaluate'))})

setMethod("evaluate", "Algorithm.SDM", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env, ...) {
  # Parameters
  text.cv.param = character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param = paste0(text.cv.param,'|',cv.param[i])
  }
  obj@parameters$cv = cv
  obj@parameters$cv.param = text.cv.param
  obj@parameters$metric = metric

  if(all(obj@data$Presence %in% c(0,1))){
    # Binary data of SDM model

    # Cross-validation
    metric = switch(metric,
                    'Kappa' = 'maxKappa',
                    'CCR' = 'max.prop.correct',
                    'TSS' = 'max.sensitivity+specificity',
                    'SES' = 'sensitivity=specificity',
                    'LW' = 'min.occurence.prediction',
                    'ROC' = 'min.ROC.plot.distance')
    if(cv == 'holdout') {
      for (i in 1:cv.param[2]) {
        data = obj@data
        data$train = FALSE
        for (p in 0:1) {
          datap = data[which(data$Presence == p),]
          datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[1]))] = TRUE
          data[which(data$Presence == p),] = datap
        }
        evaldata = data[-which(data$train),]
        evaldata = evaldata[-which(names(data) == 'train')]
        traindata = data[which(data$train),]
        traindata = traindata[-which(names(data) == 'train')]
        trainobj = obj
        trainobj@data = traindata
        model = get_model(trainobj, ...)
        predicted.values = predict(model, evaldata)
        threshold = optim.thresh(evaldata$Presence, predicted.values, thresh)
        threshold = mean(threshold[[which(names(threshold) == metric)]])
        roweval = accuracy(evaldata$Presence, predicted.values, threshold)
        if(i==1) {
          evaluation = roweval
        } else {
          evaluation = rbind(evaluation, roweval)
        }
      }
    } else if(cv == 'k-fold') {
      k = cv.param[1]
      rep = cv.param[2]
      for (i in seq_len(length(rep))) {
        data = obj@data
        data$fold = 0
        for (p in 0:1) {
          datap = data[which(data$Presence == p),]
          indices = seq_len(length(datap$fold))
          fold = 1
          while(length(indices) > 0){
            j = sample(indices, 1)
            datap$fold[j] = fold
            indices = indices[-which(indices == j)]
            if(fold != k) {fold = fold + 1} else {fold = 1}
          }
          data[which(data$Presence == p),] = datap
        }
        for(j in 1:k) {
          evaldata = data[which(data$fold == j),]
          evaldata = evaldata[-which(names(data) == 'fold')]
          traindata = data[which(data$fold != j),]
          traindata = traindata[-which(names(data) == 'fold')]
          trainobj = obj
          trainobj@data = traindata
          model = get_model(trainobj, ...)
          predicted.values = predict(model, evaldata)
          threshold = optim.thresh(evaldata$Presence, predicted.values, thresh)
          threshold = mean(threshold[[which(names(threshold) == metric)]])
          roweval = accuracy(evaldata$Presence, predicted.values, threshold)
          if(i==1 && j==1) {
            evaluation = roweval
          } else {
            evaluation = rbind(evaluation, roweval)
          }
        }
      }
    } else if(cv == 'LOO') {
      data = obj@data
      predicted.values = c()
      for(j in seq_len(length(data[,1]))) {
        evaldata = data[j,]
        traindata = data[-j,]
        trainobj = obj
        trainobj@data = traindata
        model = get_model(trainobj, ...)
        predicted.values[j] = predict(model, evaldata)
      }
      threshold = optim.thresh(data$Presence, predicted.values, thresh)
      threshold = mean(threshold[[which(names(threshold) == metric)]])
      evaluation = accuracy(data$Presence, predicted.values, threshold)
    }
    obj@evaluation = evaluation[1,]
    for(i in seq_len(length(evaluation))) {
      obj@evaluation[i] = mean(evaluation[,i], na.rm = TRUE)
    }

  } else {
    # Conitnuous values of MEMs
    warning('Evaluation is not yet implemented for continuous data of MEMs !')
  }

  return(obj)})

setMethod("evaluate", "MAXENT.SDM", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env, ...) {
  # Parameters
  text.cv.param = character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param = paste0(text.cv.param,'|',cv.param[i])
  }
  obj@parameters$cv = cv
  obj@parameters$cv.param = text.cv.param
  obj@parameters$metric = metric

  if(all(obj@data$Presence %in% c(0,1))){
    # Binary data of SDM model

    metric = switch(metric,
                    'Kappa' = 'maxKappa',
                    'CCR' = 'max.prop.correct',
                    'TSS' = 'max.sensitivity+specificity',
                    'SES' = 'sensitivity=specificity',
                    'LW' = 'min.occurence.prediction',
                    'ROC' = 'min.ROC.plot.distance')
    if(cv == 'holdout') {
      for (i in 1:cv.param[2]) {
        data = obj@data
        data$train = FALSE
        for (p in 0:1) {
          datap = data[which(data$Presence == p),]
          datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[1]))] = TRUE
          data[which(data$Presence == p),] = datap
        }
        evaldata = data[-which(data$train),]
        evaldata = evaldata[-which(names(data) == 'train')]
        traindata = data[which(data$train),]
        traindata = traindata[-which(names(data) == 'train')]
        trainobj = obj
        trainobj@data = traindata
        model = get_model(trainobj, Env)
        predicted.values = predict(model, evaldata)
        threshold = optim.thresh(evaldata$Presence, predicted.values, thresh)
        threshold = mean(threshold[[which(names(threshold) == metric)]])
        roweval = accuracy(evaldata$Presence, predicted.values, threshold)
        if(i==1) {
          evaluation = roweval
        } else {
          evaluation = rbind(evaluation, roweval)
        }
      }
    } else {
      if(cv == 'LOO') {
        k = length(obj@data$Presence)
        rep = cv.param[1]
      }
      if(cv == 'k-fold') {
        k = cv.param[1]
        rep = cv.param[2]
      }
      for (i in seq_len(length(rep))) {
        data = obj@data
        data$fold = 0
        for (p in 0:1) {
          datap = data[which(data$Presence == p),]
          indices = seq_len(length(datap$fold))
          fold = 1
          while(length(indices) > 0){
            j = sample(indices, 1)
            datap$fold[j] = fold
            indices = indices[-which(indices == j)]
            if(fold != k) {fold = fold + 1} else {fold = 1}
          }
          data[which(data$Presence == p),] = datap
        }
        for(j in 1:k) {
          evaldata = data[which(data$fold == j),]
          evaldata = evaldata[-which(names(data) == 'fold')]
          traindata = data[which(data$fold != j),]
          traindata = traindata[-which(names(data) == 'fold')]
          trainobj = obj
          trainobj@data = traindata
          model = get_model(trainobj, Env)
          predicted.values = predict(model, evaldata)
          threshold = optim.thresh(evaldata$Presence, predicted.values, thresh)
          threshold = mean(threshold[[which(names(threshold) == metric)]])
          roweval = accuracy(evaldata$Presence, predicted.values, threshold)
          if(i==1 && j==1) {
            evaluation = roweval
          } else {
            evaluation = rbind(evaluation, roweval)
          }
        }
      }
    }
    obj@evaluation = evaluation[1,]
    for(i in seq_len(length(evaluation))) {
      obj@evaluation[i] = mean(evaluation[,i], na.rm = TRUE)
    }

  }  else {
    # Conitnuous values of MEMs
    warning('Evaluation is not yet implemented for continuous data of MEMs !')
  }

  return(obj)})
