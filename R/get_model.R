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
#' @importFrom utils capture.output
NULL

setGeneric("get_model", function(obj, ...) {
  return(standardGeneric("get_model"))
})

setMethod("get_model", "Algorithm.SDM", function(obj, ....) {
  return(obj)
})

setMethod("get_model", "GLM.SDM", function(obj, glm.args=list(), ...) {
  # set defaults if arguments are NULL
  if(is.null(glm.args$test))
    glm.args$test <- "AIC"

  if(is.null(glm.args$control))
    glm.args$control <- glm.control(epsilon = 1e-08, maxit = 500)

  # information on custom settings, could be saved somewhere in future versions
  # arg.settings <- unlist(glm.args)
  # arg.settings <- paste(paste0(names(arg.settings),"=",arg.settings),collapse=",")

  if(is.null(glm.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    glm.args$data <- data
  } else {data <- glm.args$data}

  if(is.null(glm.args$formula)){
    formula <- "Presence ~"
    for (i in 2:length(data)) {
      var <- names(data[i])
      if (i != 2) {
        formula <- paste(formula, "+", var)
      } else {
        formula <- paste(formula, var)
      }
    }
    glm.args$formula <- formula(formula)
  }

  # call GLM
  model <- do.call(glm, glm.args)

  for (i in seq_len(length(data))) {
    if (is.factor(data[, i])) {
      model$xlevels[[which(names(model$xlevels) == paste0(names(data)[i]))]] <- levels(data[,
                                                                                            i])
    }
  }
  return(model)
})

setMethod("get_model", "GAM.SDM", function(obj, gam.args=list(), ...) {
  # set defaults if arguments are NULL
  if(is.null(gam.args$test))
    gam.args$test <- "AIC"

  if(is.null(gam.args$control))
    gam.args$control <- gam.control(epsilon = 1e-08,maxit = 500)

  # information on custom settings, could be saved somewhere in future versions
  # arg.settings <- unlist(gam.args)
  # arg.settings <- paste(paste0(names(arg.settings),"=",arg.settings),collapse=",")

  if(is.null(gam.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    gam.args$data <- data
  } else {data <- gam.args$data}

  if(is.null(gam.args$formula)){
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
    gam.args$formula <- formula(formula)
  }

  # call GAM
  model <- do.call(gam, gam.args)

  for (i in seq_len(length(data))) {
    if (is.factor(data[, i])) {
      model$xlevels[[which(names(model$xlevels) == names(data)[i])]] <- levels(data[,
                                                                                    i])
    }
  }
  return(model)
})

setMethod("get_model", "MARS.SDM", function(obj, mars.args=list(), ...) {
  # set defaults if arguments are NULL
  # careful, earth requires a certain order of arguments (formula,data,etc.)
  if(is.null(mars.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    mars.args$data <- data
  }
  if(is.null(mars.args$degree))
    mars.args$degree <- 2
  if(is.null(mars.args$formula)){
    mars.args <- c(formula=Presence ~ ., mars.args)
  } else {
    mars.args <- c(mars.args$formula,mars.args[-which(names(mars.args)=="formula")])
  }

  # call MARS
  model <- do.call(earth, mars.args)
  return(model)
})

setMethod("get_model", "CTA.SDM", function(obj, cta.args=list(), ...) {
  # set defaults if arguments are NULL
  if(is.null(cta.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    cta.args$data <- data
  }
  if(is.null(cta.args$control))
    cta.args$control <- rpart.control(minbucket = 1, xval = 3)
  if(is.null(cta.args$formula)){
    cta.args <- c(formula=Presence ~ ., cta.args)
    } else {
  cta.args <- c(cta.args$formula,cta.args[-which(names(cta.args)=="formula")])
  }

  # call CTA
  model <- do.call(rpart, cta.args)

  return(model)
})


setMethod("get_model", "GBM.SDM", function(obj, gbm.args=list(), ...) {

  # set defaults if arguments are NULL
  if(is.null(gbm.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    gbm.args$data <- data
  } else {data <- gbm.args$data}
  if(is.null(gbm.args$distribution)){
    if (all(data$Presence %in% c(0, 1)))
      distribution  <- "bernoulli"
    else
      distribution <- "gaussian"
    gbm.args$distribution <- distribution
  }

  if(is.null(gbm.args$n.minobsinnode))
    gbm.args$n.minobsinnode <- 1
  if(is.null(gbm.args$shrinkage))
    gbm.args$shrinkage <- 0.001
  if(is.null(gbm.args$bag.fraction))
    gbm.args$bag.fraction <- 0.5
  if(is.null(gbm.args$n.cores))
    gbm.args$n.cores <- 1
  if(is.null(gbm.args$train.fraction))
    gbm.args$train.fraction <- 1
  if(is.null(gbm.args$cv.folds))
    gbm.args$cv.folds <- 3
  if(is.null(gbm.args$n.trees))
    gbm.args$n.trees <- 2500
  if(is.null(gbm.args$verbose))
    gbm.args$verbose <- FALSE

  if(is.null(gbm.args$formula)){
    gbm.args <- c(formula=Presence ~ ., gbm.args)
  } else {
    gbm.args <- c(gbm.args$formula,gbm.args[-which(names(gbm.args)=="formula")])
  }

  # call GBM
  capture.output(model <- do.call(gbm, gbm.args), file = nullfile())
  return(model)
})

setMethod("get_model", "RF.SDM", function(obj, rf.args=list(), ...) {

  # set defaults if arguments are NULL
  if(is.null(rf.args$data)){
    data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
    rf.args$data <- data
  }

  if(is.null(rf.args$do.classif))
    rf.args$do.classif <- TRUE
  if(is.null(rf.args$ntree))
    rf.args$ntree <- 2500
  if(is.null(rf.args$nodesize))
    rf.args$nodesize <- 1

  if(is.null(rf.args$formula)){
    rf.args <- c(formula=Presence ~ ., rf.args)
  } else {
    rf.args <- c(rf.args$formula,rf.args[-which(names(rf.args)=="formula")])
  }

  # call RF
  model <- suppressWarnings(do.call(randomForest, rf.args))
  return(model)
})

setMethod("get_model", "MAXENT.SDM", function(obj, maxent.args=list(), ...) {
  # set defaults if arguments are NULL (keep order!)
  if(is.null(maxent.args$p)){
    maxent.args <- c(list(p=obj@data$Presence[which(obj@data$train)]),maxent.args)
  } else {
    maxent.args <- c(list(p=maxent.args$p), maxent.args[-which(names(maxent.args)=="p")])
  }

  if(is.null(maxent.args$x)){
    maxent.args <- c(list(x=obj@data[which(obj@data$train),which(!colnames(obj@data)%in%c("X","Y","Presence","train"))]),maxent.args)
    if(is.null(maxent.args$factors)){
      factors <- c()
      for (i in 4:length(names(obj@data))) {
        if (is.factor(obj@data[, i])) {
          factors <- c(factors, names(obj@data)[i])
        }
        maxent.args$factors <- factors
      }
    }
  } else {
    maxent.args <- c(list(x=maxent.args$x), maxent.args[-which(names(maxent.args)=="x")])
  }

  # call MAXENT
  model <- do.call(maxent, maxent.args)
  return(model)
})

setMethod("get_model", "ANN.SDM", function(obj, ann.args=list(), ...) {
  # set defaults if arguments are NULL
  if(is.null(ann.args$data)){
    ann.args$data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
  }
  if(is.null(ann.args$size))
    ann.args$size <- 6
  if(is.null(ann.args$maxit))
    ann.args$maxit <- 500
  if(is.null(ann.args$trace))
    ann.args$trace <- FALSE
  if(is.null(ann.args$formula)){
    ann.args <- c(formula=Presence ~ ., ann.args)
  } else {
    ann.args <- c(list(formula=ann.args$formula), ann.args[-which(names(ann.args)=="formula")])
  }

  # call ANN
  model <- do.call(nnet, ann.args)
  return(model)
})

setMethod("get_model", "SVM.SDM", function(obj, svm.args=list(), ...) {

  # set defaults if arguments are NULL
  if(is.null(svm.args$data)){
    svm.args$data <- obj@data[which(obj@data$train),-c(which(names(obj@data) %in% c("X","Y","train")))]
  }
  if(is.null(svm.args$type))
    svm.args$type <- "eps-regression"
  if(is.null(svm.args$gamma))
    svm.args$gamma <- 1/(length(svm.args$data) - 1)
  if(is.null(svm.args$kernel))
    svm.args$kernel <- "radial"
  if(is.null(svm.args$epsilon))
    svm.args$epsilon <- 1e-08
  if(is.null(svm.args$cross))
    svm.args$cross <- 3

  if(is.null(svm.args$formula)){
    svm.args <- c(formula=Presence ~ ., svm.args)
  } else {
    svm.args <- c(list(formula=svm.args$formula), svm.args[-which(names(svm.args)=="formula")])
  }

  # call SVM
  model <- suppressWarnings(do.call(svm, svm.args))
  return(model)
})

