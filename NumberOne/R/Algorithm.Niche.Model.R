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
#' @import nnet nnet
#' @importFrom e1071 svm
NULL

##### New generics ##### ----
setGeneric('evaluate', function(obj, thresh = 1001, metric = 'AUC') {return(standardGeneric('evaluate'))})
setGeneric('get_PA', function(obj) {return(standardGeneric('get_PA'))})
setGeneric('PA.select', function(obj, Env, ...) {return(standardGeneric('PA.select'))})
setGeneric('data.values', function(obj, Env, na.rm = T) {return(standardGeneric('data.values'))})
setGeneric('get_model', function(obj, ...) {return(standardGeneric('get_model'))})
setGeneric('project', function(obj, Env, ...) {return(standardGeneric('project'))})
setGeneric('evaluate.axes', function(obj, thresh = 1001, Env, axes.metric = 'AUC', ...) {return(standardGeneric('evaluate.axes'))})
setGeneric('ensemble', function(x, ..., name = NULL,ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = T, thresh = 1001, uncertainity = T) {return(standardGeneric('ensemble'))})
setGeneric('save.enm', function (enm, ...) {return(standardGeneric('save.enm'))})
setGeneric('save.stack', function (stack, ...) {return(standardGeneric('save.stack'))})
setGeneric('stacking', function(enm, ...) {return(standardGeneric('stacking'))})

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
setClass('Algorithm.Niche.Model',
         contains = 'Niche.Model')

# Class generator
#' @export
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

# 3 - Methods definition #
setMethod('get_PA', "Algorithm.Niche.Model", function(obj) {return(obj)})

setMethod('PA.select', "Algorithm.Niche.Model", function(obj, Env, PA = NULL, train.frac = 0.7) {
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
  data.PA$Train = F
  data.PA$Train[sample.int(length(data.PA$Presence), round(length(data.PA$Presence)*train.frac))] = T
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

setMethod('project', "Algorithm.Niche.Model",  function(obj, Env, ...) {
  model = get_model(obj, ...)
  #proj = predict(Env, model) # Previous classic function that don't take correctly factors variable into account
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

setMethod('evaluate.axes', "Algorithm.Niche.Model", function(obj, thresh = 1001, Env,
                                                             axes.metric = 'AUC', ...) {
  obj@parameters$axes.metric = axes.metric
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  if (axes.metric == 'Pearson') {
    data = obj@data[which(!obj@data$Train),]
    o.predicted.values = extract(obj@projection, data[c('X','Y')]) # original model predicted values
  }

  for (i in 5:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    data = obj.axes@data[which(!obj.axes@data$Train),]
    model = get_model(obj.axes,...)
    predicted.values = predict(model, data)
    if (axes.metric != 'Pearson') {
      threshold = optim.thresh(data$Presence, predicted.values, thresh)
      threshold = mean(threshold$`max.sensitivity+specificity`)
      evaluation = accuracy(data$Presence, predicted.values, threshold)
      obj@variables.importance[1,(i-4)] = obj@evaluation[1,which(names(obj@evaluation) == axes.metric)]  - evaluation[1,which(names(evaluation) == axes.metric)]
    } else {
      obj@variables.importance[(i-4)] = cor(predicted.values, o.predicted.values)
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
            }
            model = glm(formula(formula), data = data, test = test, control = glm.control(epsilon = epsilon, maxit = maxit))
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
              if (!is.factor(data[,i])) {formula = paste0(formula,' + s(',var,')')}
            }
            model = gam(formula(formula), data = data, test = test, control = gam.control(epsilon = epsilon, maxit = maxit))
            #             for(i in 1:length(data)) {
            #               if(is.factor(data[,i])) {
            #                 model$xlevels[[which(names(model$xlevels) == names(data)[i])]] = levels(data[,i])
            #               }
            #             }
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = rpart(Presence ~ ., data = data, method = 'class',
                          control = rpart.control(minbucket = final.leave, xval = cv))
            return(model)})

setMethod('evaluate.axes', "CTA.Niche.Model", function(obj, thresh = 1001, ...) {
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  for (i in 5:length(obj@data)) {
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    data = obj.axes@data[which(!obj.axes@data$Train),]
    model = get_model(obj.axes,...)
    predicted.values = predict(model,data)[,2]
    threshold = optim.thresh(data$Presence, predicted.values, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    evaluation = accuracy(data$Presence, predicted.values, threshold)
    obj@variables.importance[(i-4)] = obj@evaluation$AUC - evaluation$AUC
  }
  obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
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
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
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
  data = obj@data[which(obj@data$Train),]
  for(i in 4:length(names(obj@data))) {if(is.factor(obj@data[,i])) {factors = c(factors, names(obj@data)[i])}}
  model = maxent(x = Env, p = obj@data[which(data$Presence == 1),1:2],
                 a = obj@data[which(data$Presence == 0),1:2], factors = factors)
  return(model)})

setMethod('project', "MAXENT.Niche.Model", function(obj, Env, ...) {
  model = get_model(obj, Env)
  proj = predict(Env, model)
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  return(obj)})

setMethod('evaluate.axes', "MAXENT.Niche.Model", function(obj, thresh = 1001, Env, ...) {
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  for (i in 5:length(obj@data)) {
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    model = get_model(obj.axes, Env[[-(i-4)]])
    obj.axes@data = obj.axes@data[which(!obj.axes@data$Train),]
    predicted.values = predict(model, obj.axes@data)
    threshold = optim.thresh(obj.axes@data$Presence, predicted.values, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    evaluation = accuracy(obj.axes@data$Presence, predicted.values, threshold)
    obj@variables.importance[(i-4)] = obj@evaluation$AUC - evaluation$AUC
  }
  obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
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
  data = obj@data[which(obj@data$Train),]
  data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
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
  data = obj@data[which(obj@data$Train),]
  data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
  model = svm(Presence ~ ., data = data, type = 'eps-regression',
              gamma = 1/(length(data)-1), kernel = 'radial', epsilon = epsilon, cross = cv)
  return(model)})
