#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
#' @importFrom mgcv gam gam.control
#' @importFrom earth earth
#' @importFrom rpart rpart rpart.control
#' @importFrom gbm gbm
#' @importFrom randomForest randomForest
#' @importFrom dismo maxent
#' @importFrom nnet nnet
#' @importFrom e1071 svm
#' @importFrom stats aggregate.data.frame cor glm glm.control rbinom runif sd var
#' @importFrom utils lsf.str read.csv read.csv2 tail write.csv
NULL

setGeneric("get_model", function(obj, ...) {
  return(standardGeneric("get_model"))
})

setMethod("get_model", "Algorithm.SDM", function(obj, ....) {
  return(obj)
})

setMethod("get_model", "GLM.SDM", function(obj, test = "AIC", epsilon = 1e-08,
                                           maxit = 500, ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  formula <- "Presence ~"
  for (i in 2:length(data)) {
    var <- names(data[i])
    if (i != 2) {
      formula <- paste(formula, "+", var)
    } else {
      formula <- paste(formula, var)
    }
  }
  model <- glm(formula(formula), data = data, test = test, control = glm.control(epsilon = epsilon,
                                                                                 maxit = maxit))
  for (i in seq_len(length(data))) {
    if (is.factor(data[, i])) {
      model$xlevels[[which(names(model$xlevels) == paste0(names(data)[i]))]] <- levels(data[,
                                                                                            i])
    }
  }
  return(model)
})

setMethod("get_model", "GAM.SDM", function(obj, test = "AIC", epsilon = 1e-08,
                                           maxit = 500, ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  formula <- "Presence ~"
  for (i in 2:length(data)) {
    var <- names(data[i])
    if (i != 2) {
      formula <- paste(formula, "+", var)
    } else {
      formula <- paste(formula, var)
    }
    if (!is.factor(data[, i])) {
      formula <- paste0(formula, " + s(", var, ")")
    }
  }
  model <- gam(formula(formula), data = data, test = test, control = gam.control(epsilon = epsilon,
                                                                                 maxit = maxit))
  for (i in seq_len(length(data))) {
    if (is.factor(data[, i])) {
      model$xlevels[[which(names(model$xlevels) == names(data)[i])]] <- levels(data[,
                                                                                    i])
    }
  }
  return(model)
})

setMethod("get_model", "MARS.SDM", function(obj, degree = 2, ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  model <- earth(Presence ~ ., data = data, degree = 2)
  return(model)
})

setMethod("get_model", "CTA.SDM", function(obj, final.leave = 1, algocv = 3,
                                           ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  model <- rpart(Presence ~ ., data = data, control = rpart.control(minbucket = final.leave,
                                                                    xval = algocv))
  return(model)
})

setMethod("get_model", "GBM.SDM", function(obj, trees = 2500, final.leave = 1,
                                           algocv = 3, thresh.shrink = 0.001, n.cores = NULL, ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  if (all(data$Presence %in% c(0, 1)))
    distribution  <- "bernoulli"
  else
    distribution <- "gaussian"
  model <- gbm(Presence ~ ., data = data, distribution = distribution,
               n.minobsinnode = final.leave, shrinkage = thresh.shrink, bag.fraction = 0.5,
               n.cores = n.cores, train.fraction = 1, cv.folds = algocv, n.trees = trees)
  return(model)
})

setMethod("get_model", "RF.SDM", function(obj, trees = 2500, final.leave = 1,
                                          ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  model <- randomForest(Presence ~ ., data = data, do.classif = TRUE,
                        ntree = trees, nodesize = final.leave, maxnodes = NULL)
  return(model)
})

setMethod("get_model", "MAXENT.SDM", function(obj, Env, ...) {
  factors <- c()
  for (i in 4:length(names(obj@data))) {
    if (is.factor(obj@data[, i])) {
      factors <- c(factors, names(obj@data)[i])
    }
  }
  model <- maxent(x = Env, p = obj@data[which(obj@data$Presence == 1),
                                        1:2], a = obj@data[which(obj@data$Presence == 0), 1:2], factors = factors)
  return(model)
})

setMethod("get_model", "ANN.SDM", function(obj, maxit = 500, ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  model <- nnet(Presence ~ ., data = data, size = 6, maxit = maxit)
  return(model)
})

setMethod("get_model", "SVM.SDM", function(obj, epsilon = 1e-08, algocv = 3,
                                           ...) {
  data <- obj@data[-c(which(names(obj@data) == "X"), which(names(obj@data) ==
                                                             "Y"))]
  model <- svm(Presence ~ ., data = data, type = "eps-regression",
               gamma = 1/(length(data) - 1), kernel = "radial", epsilon = epsilon, cross = algocv)
  return(model)
})
