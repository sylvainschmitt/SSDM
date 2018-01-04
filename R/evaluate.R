#' @include Algorithm.SDM.R
#' @include Stacked.SDM.R
#' @import methods
#' @importFrom SDMTools optim.thresh accuracy
#' @importFrom stats aggregate.data.frame cor glm glm.control rbinom runif sd var
#' @importFrom utils lsf.str read.csv read.csv2 tail write.csv
#' @importFrom raster reclassify rasterize extract stack
#' @importFrom sp SpatialPoints coordinates
NULL

#' Evaluate
#'
#' Evaluation of SDM or ESDM habitat suitability predictions or evalaution of
#' SSDM floristic composition with Pottier et al, 2013 method (see reference
#' below)
#'
#' @param obj Stacked.SDM. SSDM to evaluate
#' @param cv character. Method of cross-validation used to evaluate the SDM (see
#'  details below).
#' @param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the SDM (see details below).
#' @param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#' @param metric character. Metric(s) used to select the best SDMs that will be
#'  included in the ensemble SDM (see details below).
#' @param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}).
#' @param ... unused argument
#'
#' @return SDM/ESDM/SSDM evaluation in a data.frame
#'
#' @examples
#'
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Evaluation
#' evaluate(SSDM)
#'
#' }
#'
#' @references Pottier, J., Dubuis, A., Pellissier, L., Maiorano, L., Rossier,
#'   L., Randin, C. F., Guisan, A. (2013). The accuracy of plant assemblage
#'   prediction from species distribution models varies along environmental
#'   gradients. Global Ecology and Biogeography, 22(1), 52-63.
#'   https://doi.org/10.1111/j.1466-8238.2012.00790.x
#'
#' @name evaluate
NULL

#' @rdname evaluate
#' @export
setGeneric('evaluate', function(obj, ...) {return(standardGeneric('evaluate'))})

#' @rdname evaluate
#' @export
setMethod("evaluate", "Algorithm.SDM", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env, ...) {
  # Parameters
  text.cv.param <- character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param <- paste0(text.cv.param, "|", cv.param[i])
  }
  obj@parameters$cv <- cv
  obj@parameters$cv.param <- text.cv.param
  obj@parameters$metric <- metric

  if (all(obj@data$Presence %in% c(0, 1))) {
    # Binary data of SDM model

    # Cross-validation
    metric <- switch(metric, Kappa = "maxKappa", CCR = "max.prop.correct",
                     TSS = "max.sensitivity+specificity", SES = "sensitivity=specificity",
                     LW = "min.occurence.prediction", ROC = "min.ROC.plot.distance")
    if (cv == "holdout") {
      for (i in 1:cv.param[2]) {
        data <- obj@data
        data$train <- FALSE
        for (p in 0:1) {
          datap <- data[which(data$Presence == p), ]
          datap$train[sample.int(length(datap$train), round(length(datap$train) *
                                                              cv.param[1]))] <- TRUE
          data[which(data$Presence == p), ] <- datap
        }
        evaldata <- data[-which(data$train), ]
        evaldata <- evaldata[-which(names(data) == "train")]
        traindata <- data[which(data$train), ]
        traindata <- traindata[-which(names(data) == "train")]
        trainobj <- obj
        trainobj@data <- traindata
        model <- get_model(trainobj, ...)
        predicted.values <- predict(model, evaldata)
        threshold <- optim.thresh(evaldata$Presence, predicted.values,
                                  thresh)
        threshold <- mean(threshold[[which(names(threshold) == metric)]])
        roweval <- accuracy(evaldata$Presence, predicted.values, threshold)
        if (i == 1) {
          evaluation <- roweval
        } else {
          evaluation <- rbind(evaluation, roweval)
        }
      }
    } else if (cv == "k-fold") {
      k <- cv.param[1]
      rep <- cv.param[2]
      for (i in seq_len(length(rep))) {
        data <- obj@data
        data$fold <- 0
        for (p in 0:1) {
          datap <- data[which(data$Presence == p), ]
          indices <- seq_len(length(datap$fold))
          fold <- 1
          while (length(indices) > 0) {
            j <- sample(indices, 1)
            datap$fold[j] <- fold
            indices <- indices[-which(indices == j)]
            if (fold != k) {
              fold <- fold + 1
            } else {
              fold <- 1
            }
          }
          data[which(data$Presence == p), ] <- datap
        }
        for (j in 1:k) {
          evaldata <- data[which(data$fold == j), ]
          evaldata <- evaldata[-which(names(data) == "fold")]
          traindata <- data[which(data$fold != j), ]
          traindata <- traindata[-which(names(data) == "fold")]
          trainobj <- obj
          trainobj@data <- traindata
          model <- get_model(trainobj, ...)
          predicted.values <- predict(model, evaldata)
          threshold <- optim.thresh(evaldata$Presence, predicted.values,
                                    thresh)
          threshold <- mean(threshold[[which(names(threshold) ==
                                               metric)]])
          roweval <- accuracy(evaldata$Presence, predicted.values,
                              threshold)
          if (i == 1 && j == 1) {
            evaluation <- roweval
          } else {
            evaluation <- rbind(evaluation, roweval)
          }
        }
      }
    } else if (cv == "LOO") {
      data <- obj@data
      predicted.values <- c()
      for (j in seq_len(length(data[, 1]))) {
        evaldata <- data[j, ]
        traindata <- data[-j, ]
        trainobj <- obj
        trainobj@data <- traindata
        model <- get_model(trainobj, ...)
        predicted.values[j] <- predict(model, evaldata)
      }
      threshold <- optim.thresh(data$Presence, predicted.values, thresh)
      threshold <- mean(threshold[[which(names(threshold) == metric)]])
      evaluation <- accuracy(data$Presence, predicted.values, threshold)
    }
    obj@evaluation <- evaluation[1, ]
    for (i in seq_len(length(evaluation))) {
      obj@evaluation[i] <- mean(evaluation[, i], na.rm = TRUE)
    }

  } else {
    # Conitnuous values of MEMs
    warning("Evaluation is not yet implemented for continuous data of MEMs !")
  }

  return(obj)
})

#' @rdname evaluate
#' @export
setMethod("evaluate", "MAXENT.SDM", function(obj, cv, cv.param, thresh = 1001, metric = 'SES', Env, ...) {
  # Parameters
  text.cv.param <- character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param <- paste0(text.cv.param, "|", cv.param[i])
  }
  obj@parameters$cv <- cv
  obj@parameters$cv.param <- text.cv.param
  obj@parameters$metric <- metric

  if (all(obj@data$Presence %in% c(0, 1))) {
    # Binary data of SDM model

    metric <- switch(metric, Kappa = "maxKappa", CCR = "max.prop.correct",
                     TSS = "max.sensitivity+specificity", SES = "sensitivity=specificity",
                     LW = "min.occurence.prediction", ROC = "min.ROC.plot.distance")
    if (cv == "holdout") {
      for (i in 1:cv.param[2]) {
        data <- obj@data
        data$train <- FALSE
        for (p in 0:1) {
          datap <- data[which(data$Presence == p), ]
          datap$train[sample.int(length(datap$train), round(length(datap$train) *
                                                              cv.param[1]))] <- TRUE
          data[which(data$Presence == p), ] <- datap
        }
        evaldata <- data[-which(data$train), ]
        evaldata <- evaldata[-which(names(data) == "train")]
        traindata <- data[which(data$train), ]
        traindata <- traindata[-which(names(data) == "train")]
        trainobj <- obj
        trainobj@data <- traindata
        model <- get_model(trainobj, Env)
        predicted.values <- predict(model, evaldata)
        threshold <- optim.thresh(evaldata$Presence, predicted.values,
                                  thresh)
        threshold <- mean(threshold[[which(names(threshold) == metric)]])
        roweval <- accuracy(evaldata$Presence, predicted.values, threshold)
        if (i == 1) {
          evaluation <- roweval
        } else {
          evaluation <- rbind(evaluation, roweval)
        }
      }
    } else {
      if (cv == "LOO") {
        k <- length(obj@data$Presence)
        rep <- cv.param[1]
      }
      if (cv == "k-fold") {
        k <- cv.param[1]
        rep <- cv.param[2]
      }
      for (i in seq_len(length(rep))) {
        data <- obj@data
        data$fold <- 0
        for (p in 0:1) {
          datap <- data[which(data$Presence == p), ]
          indices <- seq_len(length(datap$fold))
          fold <- 1
          while (length(indices) > 0) {
            j <- sample(indices, 1)
            datap$fold[j] <- fold
            indices <- indices[-which(indices == j)]
            if (fold != k) {
              fold <- fold + 1
            } else {
              fold <- 1
            }
          }
          data[which(data$Presence == p), ] <- datap
        }
        for (j in 1:k) {
          evaldata <- data[which(data$fold == j), ]
          evaldata <- evaldata[-which(names(data) == "fold")]
          traindata <- data[which(data$fold != j), ]
          traindata <- traindata[-which(names(data) == "fold")]
          trainobj <- obj
          trainobj@data <- traindata
          model <- get_model(trainobj, Env)
          predicted.values <- predict(model, evaldata)
          threshold <- optim.thresh(evaldata$Presence, predicted.values,
                                    thresh)
          threshold <- mean(threshold[[which(names(threshold) ==
                                               metric)]])
          roweval <- accuracy(evaldata$Presence, predicted.values,
                              threshold)
          if (i == 1 && j == 1) {
            evaluation <- roweval
          } else {
            evaluation <- rbind(evaluation, roweval)
          }
        }
      }
    }
    obj@evaluation <- evaluation[1, ]
    for (i in seq_len(length(evaluation))) {
      obj@evaluation[i] <- mean(evaluation[, i], na.rm = TRUE)
    }

  } else {
    # Conitnuous values of MEMs
    warning("Evaluation is not yet implemented for continuous data of MEMs !")
  }

  return(obj)})

#' @rdname evaluate
#' @export
setMethod("evaluate", "Stacked.SDM", function(obj, ...){

  # Observed composition
  obs <- lapply(obj@enms, function(x){
    x <- x@data[x@data$Presence == 1,1:2]
  })
  obs <- mapply(function(x, n){
    x['species'] <- n ; return(x)
  }, n = lapply(strsplit(names(obj@enms), '.', fixed = TRUE), function(x){x[[1]]}),
  x= obs, SIMPLIFY = FALSE)
  obs <- do.call('rbind', obs)
  obs <- with(obs, table(paste(X, Y), species))
  obs <- as.data.frame.matrix(obs)
  obs[obs > 1] <- 1
  XY <- do.call('rbind', strsplit(row.names(obs), ' '))
  row.names(obs) <- 1:dim(obs)[1]
  XY <- apply(XY, 1, as.character)
  XY <- apply(XY, 1, as.numeric)
  XY <- data.frame(XY)
  names(XY) <- c('X','Y')

  # Predicted composition
  pred <- stack(lapply(obj@enms, function(x){
    x@binary
  }))
  names(pred) <- unlist(lapply(strsplit(names(obj@enms), '.', fixed = TRUE), function(x){x[[1]]}))
  pred <- extract(pred, XY)

  # Confusion matrix
  conftab <- obs*10 + pred
  conf <- data.frame(
    TP = apply(conftab, 1, function(x){sum(length(which(x == 11)))}), # True positive
    FN = apply(conftab, 1, function(x){sum(length(which(x == 10)))}), # False negative
    FP = apply(conftab, 1, function(x){sum(length(which(x == 1)))}), # False positive
    TN = apply(conftab, 1, function(x){sum(length(which(x == 0)))}) # True negative
  )

  # Evaluation indices
  n <- dim(conf)[2] # Renmaing to follow Pottier et al, 2013
  a <- conf$TP
  b <- conf$FP
  c <- conf$FN
  d <- conf$TN
  kappa.inter <- ((a + b)*(a + c) + (c + d)*(d + b)) / n*n
  eval <- data.frame(
    'species richness error' = rowSums(pred) - rowSums(obs),
    'prediction success' = (a + d) / n,
    'kappa' = (((a + d) / n) - kappa.inter) / (1 - kappa.inter),
    'specificity' = d / (b + d),
    'sensitivity' = a / (a + c),
    'Jaccard' = a / (a + b + c)
  )

  # Result
  eval <- data.frame(
    mean = apply(eval, 2, mean, na.rm = TRUE),
    SD = apply(eval, 2, sd, na.rm = TRUE)
  )

  return(data.frame(t(eval)))
})
