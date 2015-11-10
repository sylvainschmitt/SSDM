#' @import methods
#' @importFrom sp Polygon Polygons SpatialPolygons bbox
#' @importFrom raster raster stack extract predict reclassify
#' @importFrom SDMTools optim.thresh accuracy
#' @importFrom mgcv gam gam.control
#' @importFrom earth earth
#' @importFrom rpart rpart rpart.control
#' @importFrom gbm gbm
#' @importFrom randomForest randomForest
#' @importFrom dismo maxent
#' @importFrom nnet nnet
#' @importFrom e1071 svm
NULL

##### New generics ##### ----
setGeneric('get_PA', function(obj) {return(standardGeneric('get_PA'))})
setGeneric('PA.select', function(obj, Env, ...) {return(standardGeneric('PA.select'))})
setGeneric('data.values', function(obj, Env, na.rm = T) {return(standardGeneric('data.values'))})
setGeneric('get_model', function(obj, ...) {return(standardGeneric('get_model'))})
setGeneric('evaluate', function(obj, cv = 'holdout', cv.param = c(0.7, 2), thresh = 1001, metric = 'SES', Env, ...) {return(standardGeneric('evaluate'))})
setGeneric('project', function(obj, Env, ...) {return(standardGeneric('project'))})
setGeneric('evaluate.axes', function(obj, cv = 'holdout', cv.param = c(0.7, 2), thresh = 1001,
                                     metric = 'SES', axes.metric = 'Pearson', Env, ...) {return(standardGeneric('evaluate.axes'))})

##### Niche Model Class ##### -----
# 1 - Class definition #
setClass('Niche.Model',
         representation(name = 'character',
                        projection = 'Raster',
                        evaluation = 'data.frame',
                        variables.importance = 'data.frame',
                        data = 'data.frame',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   projection = raster(),
                   evaluation = data.frame(),
                   variables.importance = data.frame(),
                   data = data.frame(),
                   parameters = data.frame()))

# 2 - Methods definition #
setMethod('print', 'Niche.Model', function(x, ...) {
  cat('Object of class :', class(x)[1],'\n')
  cat('Name :', x@name, '\n')
  cat('Projections : ',names(x@projection),'\n')
  print(x@evaluation)
  print(x@variables.importance)
  if(inherits(x, 'Ensemble.Niche.Model')) {
    cat('Uncertinity map :', names(x@uncertainity),'\n')
    print(x@algorithm.evaluation)
    print(x@algorithm.correlation)
  }
})

##### Algorithm Niche Model Class ##### -----
#'An S4 class to represent a specie distribution model of one algorithm
#'
#'This is an S4 class to represent a specie distribution model of one algorithm
#'(among generalized linear model, general additive model, multivariate
#'adpatative splines, generalized boosted regression models, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines). It can be obtain with \code{\link{Modelling}}.
#'
#'@slot name character. Name of the model (by default
#'  Specie.Algorithm.Niche.Model)
#'@slot projection raster. Habitat suitability map of the model
#'@slot evaluation data frame. Evaluation of the model (threshold, AUC, omission
#'  rate, sensitivity, specificity, correct proportion and Kappa)
#'@slot variables.importance data frame. Relative percentage of importance for each variable used in the model
#'@slot data data frame. Data used to realized the model
#'@slot parameters data frame. Parameters used to realized the model
#'
#'@seealso \linkS4class{Ensemble.Niche.Model} an S4 class for ensemble models,
#'  and \linkS4class{Stack.Species.Ensemble.Niche.Model} an S4 class for stack
#'  species enemble models.
#'
#' @export
setClass('Algorithm.Niche.Model',
         contains = 'Niche.Model')

# Class generator
Algorithm.Niche.Model <- function(algorithm = 'Algorithm',
                                  name = character(),
                                  projection = raster(),
                                  evaluation = data.frame(),
                                  variables.importance = data.frame(),
                                  data = data.frame(),
                                  parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  object.class = paste0(algorithm,'.Niche.Model')
  return(new(object.class, name = name, projection = projection, evaluation = evaluation, variables.importance = variables.importance, data = data, parameters = parameters))
}

# 3 - Internal undocumented methods #
setMethod('get_PA', "Algorithm.Niche.Model", function(obj) {return(obj)})

setMethod('PA.select', "Algorithm.Niche.Model", function(obj, Env, PA = NULL) {
  if (is.null(PA)) {
    PA = get_PA(obj)
    obj@parameters$PA = 'default'
  } else {
    obj@parameters$PA = paste0(as.character(PA$nb),'.',as.character(PA$strat))
  }

  # Mask defining
  if (PA$strat == '2nd') {
    cat('   second far selection \n')
    circles = list()
    for (i in 1:length(obj@data$X)) {
      x = obj@data$X[i]
      y = obj@data$Y[i]
      pts = seq(0, 2 * pi, length.out = 100)
      xy = cbind(x + 2/60 * sin(pts), y + 2/60 * cos(pts))
      circle = Polygon(xy)
      circles[i] = circle
    }
    sc= SpatialPolygons(list(Polygons(circles, 'Circles')))
    Mask = mask(Env[[1]], sc)
  } else {
    cat('   random selection \n')
    Mask = Env[[1]]}

  # Pseudo-Absences selection
  data.PA = data.frame(matrix(nrow = 0, ncol = 2))
  names(data.PA) = c('X','Y')
  if(PA$nb < 100) {nb = PA$nb*PA$nb} else {nb = 1000}
  while (length(data.PA[,1]) < PA$nb) {
    X = runif(nb, min = bbox(Mask)[1,1], max = bbox(Mask)[1,2])
    Y = runif(nb, min = bbox(Mask)[2,1],max = bbox(Mask)[2,2])
    points = data.frame(X = X, Y = Y)
    check = extract(Mask, points)
    points = points[-which(is.na(check)),]
    data.PA = rbind(data.PA, points)
  }
  data.PA = data.PA[1:PA$nb,]
  data.PA$Presence = 0
  obj@data = rbind(obj@data, data.PA)

  return(obj)})

setMethod('data.values', "Algorithm.Niche.Model", function(obj, Env, na.rm = T) {
  values = data.frame(extract(Env, cbind(obj@data$X, obj@data$Y)))

  # Categorical variables as factor
  for (i in 1:length(Env@layers)) {
    if(Env[[i]]@data@isfactor) {
      col = which(names(values) == Env[[i]]@data@names)
      values[,col] = as.factor(values[,col])
      levels(values[,col]) = Env[[i]]@data@attributes[[1]]$ID
      if(length(Env[[i]]@data@attributes[[1]]$ID) > 100) {
        warning(paste(names(Env[[i]]), 'as more than 100 levels (', length(Env[[i]]@data@attributes[[1]]$ID), ') are you sure to consider it as a factor ?'))
        }
    }
  }

  # Tables binding
  obj@data = cbind(obj@data, values)

  # NAs removing
  if(na.rm) {
    for (i in 1:length(Env@layers)) {
      if(length(which(is.na(obj@data[i+3]))) > 0) {obj@data = obj@data[-c(which(is.na(obj@data[i+3]))),]}
    }
  }

  return(obj)})

setMethod('get_model', "Algorithm.Niche.Model", function(obj, ....) {return(obj)})

setMethod("evaluate", "Algorithm.Niche.Model", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env, ...) {
  # Parameters
  text.cv.param = character()
  for (i in 1:length(cv.param)) {
    text.cv.param = paste0(text.cv.param,'|',cv.param[i])
  }
  obj@parameters$cv = cv
  obj@parameters$cv.param = text.cv.param
  obj@parameters$metric = metric
  # Cross validation
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
      data$train = F
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[1]))] = T
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
  } else {
    if(cv == 'LOO') {
      k = length(obj@data$Presence)
      rep = cv.param[1]
    }
    if(cv == 'k-fold') {
      k = cv.param[1]
      rep = cv.param[2]
    }
    for (i in 1:length(rep)) {
      data = obj@data
      data$fold = 0
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        indices = c(1:length(datap$fold))
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
  }
  obj@evaluation = evaluation[1,]
  for(i in 1:length(evaluation)){
    obj@evaluation[i] = mean(evaluation[,i], na.rm = T)
  }
  return(obj)})

setMethod('project', "Algorithm.Niche.Model",  function(obj, Env, ...) {
  model = get_model(obj, ...)
  proj = raster::predict(Env, model,
                         fun = function(model, x){
                           x= as.data.frame(x)
                           for (i in 1:length(Env@layers)) {
                             if(Env[[i]]@data@isfactor) {
                               x[,i] = as.factor(x[,i])
                               x[,i] = droplevels(x[,i])
                               levels(x[,i]) = Env[[i]]@data@attributes[[1]]$ID
                             }
                           }
                           return(predict(model, x))
                         })
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  return(obj)})

setMethod('evaluate.axes', "Algorithm.Niche.Model", function(obj, cv, cv.param, thresh = 1001,
                                                             metric = 'SES', axes.metric = 'Pearson', Env, ...) {
  obj@parameters$axes.metric = axes.metric
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-3)))
  names(obj@variables.importance) = names(obj@data)[4:length(obj@data)]
  if (axes.metric == 'Pearson') {
    o.predicted.values = predict(get_model(obj, ...), obj@data) # original model predicted values
  }

  for (i in 4:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    if (axes.metric != 'Pearson') {
      obj.axes = evaluate(obj.axes, cv, cv.param, thresh, metric)
      obj@variables.importance[1,(i-3)] = obj@evaluation[1,which(names(obj@evaluation) == axes.metric)]  - obj.axes@evaluation[1,which(names(obj.axes@evaluation) == axes.metric)]
    } else {
      predicted.values = predict(get_model(obj.axes, ...), obj.axes@data)
      obj@variables.importance[(i-3)] = cor(predicted.values, o.predicted.values)
    }
  }

  # Variable importance normalization (%)
  if(sum(obj@variables.importance) == 0) {
    all.null = T
    for (i in 1:length(obj@variables.importance[1,])) {if(obj@variables.importance[1,i] != 0) {all.null = F}}
    if (all.null) {
      obj@variables.importance[1,] = 100 / length(obj@variables.importance)
    } else {
      obj@variables.importance = obj@variables.importance * 100
    }
  } else {
    obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  }
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

##### GLM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GLM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GLM.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = 1000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GLM.Niche.Model",
          function(obj, test = 'AIC', epsilon = 1e-08, maxit = 500) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
            }
            model = glm(formula(formula), data = data, test = test,
                        control = glm.control(epsilon = epsilon, maxit = maxit))
            for(i in 1:length(data)) {
              if(is.factor(data[,i])) {
                model$xlevels[[which(names(model$xlevels) == paste0(names(data)[i]))]] = levels(data[,i])
              }
            }
            return(model)})

##### GAM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GAM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GAM.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = 1000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GAM.Niche.Model",
          function(obj, test = 'AIC', epsilon = 1e-08, maxit = 500) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
              if (!is.factor(data[,i])) {formula = paste0(formula,' + s(',var,')')}
            }
            model = gam(formula(formula), data = data, test = test,
                        control = gam.control(epsilon = epsilon, maxit = maxit))
                        for(i in 1:length(data)) {
                          if(is.factor(data[,i])) {
                            model$xlevels[[which(names(model$xlevels) == names(data)[i])]] = levels(data[,i])
                          }
                        }
            return(model)})

##### MARS Niche Model Class ##### -----

# 1 - Class definition #
setClass('MARS.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "MARS.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = 100
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "MARS.Niche.Model",
          function(obj, degree = 2) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            model = earth(Presence ~ ., data = data, degree = 2)
            return(model)})

##### CTA Niche Model Class ##### -----

# 1 - Class definition #
setClass('CTA.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "CTA.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "CTA.Niche.Model",
          function(obj, final.leave = 1, cv = 3) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            model = rpart(Presence ~ ., data = data, method = 'class',
                          control = rpart.control(minbucket = final.leave, xval = cv))
            return(model)})

setMethod("evaluate", "CTA.Niche.Model", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', ...) {
  # Need to redefine method because CTA predict method doesn't return the same thing as other algorithms predict method !
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
      data$train = F
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[1]))] = T
        data[which(data$Presence == p),] = datap
      }
      evaldata = data[-which(data$train),]
      evaldata = evaldata[-which(names(data) == 'train')]
      traindata = data[which(data$train),]
      traindata = traindata[-which(names(data) == 'train')]
      trainobj = obj
      trainobj@data = traindata
      model = get_model(trainobj, ...)
      predicted.values = predict(model, evaldata)[,1]
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
    for (i in 1:length(rep)) {
      data = obj@data
      data$fold = 0
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        indices = c(1:length(datap$fold))
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
        predicted.values = predict(model, evaldata)[,1]
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
  for(i in 1:length(evaluation)){
    obj@evaluation[i] = mean(evaluation[,i], na.rm = T)
  }
  return(obj)})

setMethod('evaluate.axes', "CTA.Niche.Model", function(obj, cv, cv.param, thresh = 1001,
                                                             metric = 'SES', axes.metric = 'Pearson', ...) {
  # Need to redefine method because CTA predict method doesn't return the same thing as other algorithms predict method !
  obj@parameters$axes.metric = axes.metric
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-3)))
  names(obj@variables.importance) = names(obj@data)[4:length(obj@data)]
  if (axes.metric == 'Pearson') {
    o.predicted.values = predict(get_model(obj, ...), obj@data)[1,] # original model predicted values
  }

  for (i in 4:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    if (axes.metric != 'Pearson') {
      obj.axes = evaluate(obj.axes, cv, cv.param, thresh, metric)
      obj@variables.importance[1,(i-3)] = obj@evaluation[1,which(names(obj@evaluation) == axes.metric)]  - obj.axes@evaluation[1,which(names(obj.axes@evaluation) == axes.metric)]
    } else {
      predicted.values = predict(get_model(obj.axes, ...), obj.axes@data)[1,]
      obj@variables.importance[(i-3)] = cor(predicted.values, o.predicted.values)
    }
  }

  # Variable importance normalization (%)
  if(sum(obj@variables.importance) == 0) {
    all.null = T
    for (i in 1:length(obj@variables.importance[1,])) {if(obj@variables.importance[1,i] != 0) {all.null = F}}
    if (all.null) {
      obj@variables.importance[1,] = 100 / length(obj@variables.importance)
    } else {
      obj@variables.importance = obj@variables.importance * 100
    }
  } else {
    obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  }
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

##### GBM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GBM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GBM.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GBM.Niche.Model",
          function(obj, trees = 2500, final.leave = 1, cv = 3, thresh.shrink = 1e-03) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            model = gbm(Presence ~ ., data = data,
                        distribution = 'bernoulli', n.minobsinnode = final.leave,
                        shrinkage = thresh.shrink, bag.fraction = 0.5,
                        train.fraction = 1, cv.folds = cv, n.trees = trees)
            return(model)})

##### RF Niche Model Class ##### -----

# 1 - Class definition #
setClass('RF.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "RF.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "RF.Niche.Model",
          function(obj, trees = 2500, final.leave = 1) {
            data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
            model = randomForest(Presence ~ ., data = data, do.classif = TRUE, ntree = trees,
                                 nodesize = final.leave, maxnodes = NULL)
            return(model)})

##### MAXENT Niche Model Class ##### -----

# 1 - Class definition #
setClass('MAXENT.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "MAXENT.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = 10000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "MAXENT.Niche.Model", function(obj, Env) {
  factors = c()
  for(i in 4:length(names(obj@data))) {if(is.factor(obj@data[,i])) {factors = c(factors, names(obj@data)[i])}}
  model = maxent(x = Env, p = obj@data[which(obj@data$Presence == 1),1:2],
                 a = obj@data[which(obj@data$Presence == 0),1:2], factors = factors)
  return(model)})

setMethod("evaluate", "MAXENT.Niche.Model", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env) {
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
      data$train = F
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[1]))] = T
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
    for (i in 1:length(rep)) {
      data = obj@data
      data$fold = 0
      for (p in 0:1) {
        datap = data[which(data$Presence == p),]
        indices = c(1:length(datap$fold))
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
  for(i in 1:length(evaluation)){
    obj@evaluation[i] = mean(evaluation[,i], na.rm = T)
  }
  return(obj)})

setMethod('project', "MAXENT.Niche.Model",  function(obj, Env) {
  model = get_model(obj, Env)
  proj = raster::predict(Env, model,
                         fun = function(model, x){
                           x= as.data.frame(x)
                           for (i in 1:length(Env@layers)) {
                             if(Env[[i]]@data@isfactor) {
                               x[,i] = as.factor(x[,i])
                               x[,i] = droplevels(x[,i])
                               levels(x[,i]) = Env[[i]]@data@attributes[[1]]$ID
                             }
                           }
                           return(predict(model, x))
                         })
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  return(obj)})

setMethod('evaluate.axes', "MAXENT.Niche.Model", function(obj, cv, cv.param, thresh = 1001,
                                                             metric = 'SES', axes.metric = 'Pearson', Env) {
  obj@parameters$axes.metric = axes.metric
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-3)))
  names(obj@variables.importance) = names(obj@data)[4:length(obj@data)]
  if (axes.metric == 'Pearson') {
    o.predicted.values = predict(get_model(obj, Env), obj@data) # original model predicted values
  }

  for (i in 4:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    if (axes.metric != 'Pearson') {
      obj.axes = evaluate(obj.axes, cv, cv.param, thresh, metric, Env[[-(i-3)]])
      obj@variables.importance[1,(i-3)] = obj@evaluation[1,which(names(obj@evaluation) == axes.metric)]  - obj.axes@evaluation[1,which(names(obj.axes@evaluation) == axes.metric)]
    } else {
      predicted.values = predict(get_model(obj.axes, Env[[-(i-3)]]), obj.axes@data)
      obj@variables.importance[(i-3)] = cor(predicted.values, o.predicted.values)
    }
  }

  # Variable importance normalization (%)
  if(sum(obj@variables.importance) == 0) {
    all.null = T
    for (i in 1:length(obj@variables.importance[1,])) {if(obj@variables.importance[1,i] != 0) {all.null = F}}
    if (all.null) {
      obj@variables.importance[1,] = 100 / length(obj@variables.importance)
    } else {
      obj@variables.importance = obj@variables.importance * 100
    }
  } else {
    obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  }
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

##### ANN Niche Model Class ##### -----

# 1 - Class definition #
setClass('ANN.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "ANN.Niche.Model",
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "ANN.Niche.Model", function(obj, maxit = 500) {
  data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
  model = nnet(Presence ~ ., data = data, size = 6, maxit = maxit)
  return(model)})

##### SVM Niche Model Class ##### -----

# 1 - Class definition #
setClass('SVM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "SVM.Niche.Model", function(obj) {
  PA = list()
  PA['nb'] = length(obj@data$Presence)
  PA['strat'] = 'random'
  return(PA)})

setMethod('get_model', "SVM.Niche.Model", function(obj, epsilon = 1e-08, cv = 3) {
  data = obj@data[-c(which(names(obj@data) == 'X'),which(names(obj@data) == 'Y'))]
  model = svm(Presence ~ ., data = data, type = 'eps-regression',
              gamma = 1/(length(data)-1), kernel = 'radial', epsilon = epsilon, cross = cv)
  return(model)})
