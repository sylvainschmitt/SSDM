#' @include Algorithm.SDM.R
#' @include Stacked.SDM.R
#' @include optim.thresh.R
#' @include accuracy.R
#' @import methods
#' @importFrom dismo evaluate threshold
#' @importFrom stats aggregate.data.frame cor glm glm.control rbinom runif sd var
#' @importFrom utils lsf.str read.csv read.csv2 tail write.csv
#' @importFrom raster reclassify rasterize extract stack
#' @importFrom sdm calibration
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
#' @param final.fit.data strategy used for fitting the final model to be returned: 'holdout'= use same train and test data as in (last) evaluation, 'all'= train model with all data (i.e. no test data) or numeric (0-1)= sample a custom training fraction (left out fraction is set aside as test data)
#' @param bin.thresh character. Classification threshold (\code{\link[dismo]{threshold}}) used to binarize model predictions into presence/absence and compute the confusion matrix (including related scores such as TPR, TNR, omission rate, Kappa, etc.).
#' @param metric (deprecated) character. Classification threshold (\code{SDMTools::optim.thresh}) used to binarize model predictions into presence/absence and compute the confusion matrix (including related scores such as TPR, TNR, omission rate, Kappa, etc.).
#' @param thresh (deprecated) integer. Number of equally spaced thresholds in the interval 0-1 (\code{SDMTools::optim.thresh}).
#' @param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}).
#' @param ... arguments for internal use (get_model), such as argument lists to be passed to the source functions (e.g. glm.args=list(test="AIC",singular.ok=FALSE))
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
#'   L., Randin, C. F., Guisan, A. (2013). The .accuracy of plant assemblage
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
setMethod("evaluate", "Algorithm.SDM", function(obj, cv, cv.param, final.fit.data='all', bin.thresh = 'SES', metric = NULL, thresh=1001, Env, ...) {
  # Parameters
  text.cv.param <- character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param <- paste0(text.cv.param, "|", cv.param[i])
  }
  obj@parameters$cv <- cv
  obj@parameters$cv.param <- text.cv.param
  if(!is.null(metric)){
    warning("Argument 'metric' is deprecated and will be removed in future versions. Please consider using 'bin.thresh' instead.")
    obj@parameters$metric <- metric
  } else {
    obj@parameters$metric <- bin.thresh
  }


  if (all(obj@data$Presence %in% c(0, 1))) {
    # Binary data of SDM model

    # translate thresholding terms (metric only for backwards compatibility)
    if(!is.null(metric)){
      metric <- switch(metric, Kappa = "maxKappa", CCR = "max.prop.correct",
                       TSS = "max.sensitivity+specificity", SES = "sensitivity=specificity",
                       LW = "min.occurence.prediction", ROC = "min.ROC.plot.distance")
    } else {
      bin.thresh <- switch(bin.thresh, Kappa = "kappa", NOM = "no_omission",
                           TSS = "spec_sens", SES = "equal_sens_spec",
                           EP = "prevalence")
    }


    # Cross-validation

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
        eval.testdata <- data[-which(data$train), ]
        eval.testdata <- eval.testdata[-which(names(data) == "train")]
        eval.traindata <- data[which(data$train), ]
        # eval.traindata <- eval.traindata[-which(names(data) == "train")]
        evalobj <- obj
        evalobj@data <- eval.traindata
        model <- get_model(evalobj, ...)
        predicted.values <- c(predict(model, eval.testdata))
        if(!is.null(metric)){
          threshval <- .optim.thresh(eval.testdata$Presence, predicted.values, thresh)
          threshval <- mean(threshval[[which(names(threshval) == metric)]])
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)],tr= threshval)
        } else {
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)])
          threshval <- dismo::threshold(roweval,stat=bin.thresh)
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], tr=threshval)
        }
        caleval <- sdm::calibration(eval.testdata$Presence,predicted.values, nbin=20, weight=TRUE)

        evaldf <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
        if (i == 1) {
          evaluation <- evaldf
        } else {
          evaluation <- rbind(evaluation, evaldf)
        }
      }
    } else if (cv == "k-fold") {
      if(cv.param[1]<2){
        warning("less than 2 folds selected, automatic adjustment to k=5")
        cv.param[1] <- 5
      }
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
            j <- ifelse(length(indices)==1,indices,sample(indices, 1))
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
          eval.testdata <- data[which(data$fold == j), ]
          eval.testdata <- eval.testdata[-which(names(data) == "fold")]
          eval.traindata <- data[which(data$fold != j), ]
          eval.traindata <- eval.traindata[-which(names(data) == "fold")]
          eval.traindata$train <- TRUE
          evalobj <- obj
          evalobj@data <- eval.traindata
          model <- get_model(evalobj, ...)
          predicted.values <- c(predict(model, eval.testdata))
          if(!is.null(metric)){
            threshval <- .optim.thresh(eval.testdata$Presence, predicted.values, thresh)
            threshval <- mean(threshval[[which(names(threshval) == metric)]])
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
          } else {
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)])
            threshval <- dismo::threshold(roweval,stat=bin.thresh)
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
          }
          caleval <- sdm::calibration(eval.testdata$Presence,predicted.values, nbin=20, weight=TRUE)

          evaldf <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
          if (i == 1 && j == 1) {
            evaluation <- evaldf
          } else {
            evaluation <- rbind(evaluation, evaldf)
          }
        }
      }
    } else if (cv == "LOO") {
      data <- obj@data
      predicted.values <- c()
      for (j in seq_len(length(data[, 1]))) {
        eval.testdata <- data[j, ]
        eval.traindata <- data[-j, ]
        eval.traindata$train <- TRUE
        evalobj <- obj
        evalobj@data <- eval.traindata
        model <- get_model(evalobj, ...)
        predicted.values[j] <- c(predict(model, eval.testdata))
      }
      if(!is.null(metric)){
        threshval <- .optim.thresh(data$Presence, predicted.values, thresh)
        threshval <- mean(threshval[[which(names(threshval) == metric)]])
        roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)], threshval)
      } else {
        roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)])
        threshval <- dismo::threshold(roweval,stat=bin.thresh)
        roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)], threshval)
      }
      caleval <- sdm::calibration(data$Presence,predicted.values, nbin=20, weight=TRUE)

      evaluation <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
    }
    obj@evaluation <- evaluation[1, ]
    for (i in seq_len(length(evaluation))) {
      obj@evaluation[i] <- mean(evaluation[, i], na.rm = TRUE)
    }

  } else {
    # Continuous values of MEMs
    warning("Evaluation is not yet implemented for continuous data of MEMs !")
  }

  # assign train/test fractions for final model training

  if(final.fit.data=='holdout'){
    obj@data <- data
  }
  if(is.numeric(final.fit.data)){
    if(final.fit.data>0 & final.fit.data<=1){
      data <- obj@data
      data$train <- FALSE
      for (p in 0:1) {
        datap <- data[which(data$Presence == p), ]
        datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[3]))] <- TRUE
        data[which(data$Presence == p), ] <- datap
      }
      obj@data <- data
    } else{
      warning("Training fraction needs to be between 0 and 1 for sampling, assuming 1")
      final.fit.data <- 'all'}
  }
  if(final.fit.data=='all'){
    data$train <- TRUE
    obj@data <- data
  }


  return(obj)
})

#' @rdname evaluate
#' @export
setMethod("evaluate", "MAXENT.SDM", function(obj, cv, cv.param, final.fit.data='all', bin.thresh = 'SES', metric = NULL, thresh = 1001, Env, ...) {
  # Parameters
  text.cv.param <- character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param <- paste0(text.cv.param, "|", cv.param[i])
  }
  obj@parameters$cv <- cv
  obj@parameters$cv.param <- text.cv.param

  if(!is.null(metric)){
    warning("Argument 'metric' is deprecated and will be removed in future versions. Please consider using 'bin.thresh' instead.")
    obj@parameters$metric <- metric
  } else {
    obj@parameters$metric <- bin.thresh
  }

  if (all(obj@data$Presence %in% c(0, 1))) {
    # Binary data of SDM model

    # translate thresholding terms (metric only for backwards compatibility)
    if(!is.null(metric)){
      metric <- switch(metric, Kappa = "maxKappa", CCR = "max.prop.correct",
                       TSS = "max.sensitivity+specificity", SES = "sensitivity=specificity",
                       LW = "min.occurence.prediction", ROC = "min.ROC.plot.distance")
    } else {
      bin.thresh <- switch(bin.thresh, Kappa = "kappa", CCR = "no_omission",
                           TSS = "spec_sens", SES = "equal_sens_spec",
                           EP = "prevalence")
    }

    if (cv == "holdout") {
      for (i in 1:cv.param[2]) {
        data <- obj@data
        data$train <- FALSE
        # only sample test presences, not test background points
        for (p in 0:1) {
          datap <- data[which(data$Presence == p), ]
          if(p==0){
            datap$train <- TRUE
            bgdata <- datap
          } else {
            datap$train[sample.int(length(datap$train), round(length(datap$train) * cv.param[1]))] <- TRUE
          }
          data[which(data$Presence == p), ] <- datap
        }
        eval.testdata <- data[-which(data$train), ]
        # add background points as test absences (can be the same as training, if not spatially biased)
        eval.testdata <- rbind(eval.testdata,bgdata)
        eval.testdata <- eval.testdata[-which(names(data) == "train")]
        eval.traindata <- data[which(data$train), ]

        evalobj <- obj
        evalobj@data <- eval.traindata
        model <- get_model(evalobj, ...)
        predicted.values <- c(predict(model, eval.testdata))
        if(!is.null(metric)){
          threshval <- .optim.thresh(eval.testdata$Presence, predicted.values, thresh)
          threshval <- mean(threshval[[which(names(threshval) == metric)]])
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
        } else {
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)])
          threshval <- dismo::threshold(roweval,stat=bin.thresh)
          roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
        }
        caleval <- sdm::calibration(eval.testdata$Presence,predicted.values, nbin=20, weight=TRUE)

        evaldf <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
        if (i == 1) {
          evaluation <- evaldf
        } else {
          evaluation <- rbind(evaluation, evaldf)
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
          if(p==0){
            bgdata <- data[which(data$Presence == p), ]
          } else {
            datap <- data[which(data$Presence == p), ]
            indices <- seq_len(length(datap$fold))
            fold <- 1
            while (length(indices) > 0) {
              j <- ifelse(length(indices)==1,indices,sample(indices, 1))
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
        }
        for (j in 1:k) {
          eval.testdata <- data[which(data$fold == j), ]
          eval.testdata <- rbind(eval.testdata,bgdata)
          eval.testdata <- eval.testdata[-which(names(data) == "fold")]
          eval.traindata <- data[which(data$fold != j), ]
          eval.traindata <- eval.traindata[-which(names(data) == "fold")]
          eval.traindata$train <- TRUE
          evalobj <- obj
          evalobj@data <- eval.traindata
          model <- get_model(evalobj, Env)
          predicted.values <- c(predict(model, eval.testdata))
          if(!is.null(metric)){
            threshval <- .optim.thresh(eval.testdata$Presence, predicted.values, thresh)
            threshval <- mean(threshval[[which(names(threshval) == metric)]])
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], tr=threshval)
          } else {
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)])
            threshval <- dismo::threshold(roweval,stat=bin.thresh)
            roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], tr=threshval)
          }
          caleval <- sdm::calibration(eval.testdata$Presence,predicted.values, nbin=20, weight=TRUE)

          evaldf <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
          if (i == 1 && j == 1) {
            evaluation <- evaldf
          } else {
            evaluation <- rbind(evaluation, evaldf)
          }
        }
      }
    }
    obj@evaluation <- evaluation[1, ]
    for (i in seq_len(length(evaluation))) {
      obj@evaluation[i] <- mean(evaluation[, i], na.rm = TRUE)
    }

  } else {
    # Continuous values of MEMs
    warning("Evaluation is not yet implemented for continuous data of MEMs !")
  }

  # assign train/test fractions for final model training

  if(final.fit.data=='holdout'){
    obj@data <- data
  }
  if(is.numeric(final.fit.data)){
    if(final.fit.data>0 & final.fit.data<=1){
      data <- obj@data
      data$train <- FALSE
      for (p in 0:1) {
        datap <- data[which(data$Presence == p), ]
        datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[3]))] <- TRUE
        data[which(data$Presence == p), ] <- datap
      }
      obj@data <- data
    } else{
      warning("Training fraction needs to be between 0 and 1 for sampling, assuming 1")
      final.fit.data <- 'all'}
  }
  if(final.fit.data=='all'){
    data$train <- TRUE
    obj@data <- data
  }


  return(obj)})

#' @rdname evaluate
#' @export
setMethod("evaluate", "Stacked.SDM", function(obj, ...){

  # Observed composition
  obs <- lapply(obj@esdms, function(x){
    x <- x@data[x@data$Presence == 1,1:2]
  })
  obs <- mapply(function(x, n){
    x['species'] <- n ; return(x)
  }, n = lapply(strsplit(names(obj@esdms), '.', fixed = TRUE), function(x){x[[1]]}),
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
  pred <- stack(lapply(obj@esdms, function(x){
    x@binary
  }))
  names(pred) <- unlist(lapply(strsplit(names(obj@esdms), '.', fixed = TRUE), function(x){x[[1]]}))
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
  n <- dim(conf)[2] # Renaming according to Pottier et al, 2013
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
