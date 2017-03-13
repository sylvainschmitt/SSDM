#' @include Algorithm.SDM.R
#' @import methods
NULL
setGeneric('evaluate.axes', function(obj, cv = 'holdout', cv.param = c(0.7, 2), thresh = 1001,
                                     metric = 'SES', axes.metric = 'Pearson', Env, ...) {return(standardGeneric('evaluate.axes'))})

setMethod('evaluate.axes', "Algorithm.SDM", function(obj, cv, cv.param, thresh = 1001,
                                                     metric = 'SES', axes.metric = 'Pearson', Env, ...) {
  obj@parameters$axes.metric <- axes.metric
  obj@variable.importance <- data.frame(matrix(nrow = 1, ncol = (length(obj@data) -
                                                                   3)))
  names(obj@variable.importance) <- names(obj@data)[4:length(obj@data)]
  if (axes.metric == "Pearson") {
    o.predicted.values <- predict(get_model(obj, ...), obj@data)  # original model predicted values
  }
  for (i in 4:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes <- obj
    obj.axes@data <- obj.axes@data[-i]
    if (axes.metric != "Pearson") {
      obj.axes <- evaluate(obj.axes, cv, cv.param, thresh, metric)
      obj@variable.importance[1, (i - 3)] <- obj@evaluation[1, which(names(obj@evaluation) ==
                                                                       axes.metric)] - obj.axes@evaluation[1, which(names(obj.axes@evaluation) ==
                                                                                                                      axes.metric)]
    } else {
      model.axes <- get_model(obj.axes, ...)
      predicted.values <- predict(model.axes, obj.axes@data)
      c <- cor(predicted.values, o.predicted.values)
      if (is.na(c) || !is.numeric(c)) {
        c <- 0
      }
      obj@variable.importance[(i - 3)] <- 1 - c
    }
  }

  # Variable importance normalization (%)
  if (sum(obj@variable.importance) == 0) {
    all.null <- TRUE
    for (i in seq_len(length(obj@variable.importance[1, ]))) {
      if (obj@variable.importance[1, i] != 0) {
        all.null <- FALSE
      }
    }
    if (all.null) {
      obj@variable.importance[1, ] <- 100/length(obj@variable.importance)
    } else {
      obj@variable.importance <- obj@variable.importance * 100
    }
  } else {
    obj@variable.importance <- obj@variable.importance/sum(obj@variable.importance) *
      100
  }
  row.names(obj@variable.importance) <- "Axes.evaluation"
  return(obj)
})

setMethod('evaluate.axes', "MAXENT.SDM", function(obj, cv, cv.param, thresh = 1001,
                                                  metric = 'SES', axes.metric = 'Pearson', Env, ...) {
  obj@parameters$axes.metric <- axes.metric
  obj@variable.importance <- data.frame(matrix(nrow = 1, ncol = (length(obj@data) -
                                                                   3)))
  names(obj@variable.importance) <- names(obj@data)[4:length(obj@data)]
  if (axes.metric == "Pearson") {
    o.predicted.values <- predict(get_model(obj, Env), obj@data)  # original model predicted values
  }

  for (i in 4:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes <- obj
    obj.axes@data <- obj.axes@data[-i]
    if (axes.metric != "Pearson") {
      obj.axes <- evaluate(obj.axes, cv, cv.param, thresh, metric, Env[[-(i -
                                                                            3)]], ...)
      obj@variable.importance[1, (i - 3)] <- obj@evaluation[1, which(names(obj@evaluation) ==
                                                                       axes.metric)] - obj.axes@evaluation[1, which(names(obj.axes@evaluation) ==
                                                                                                                      axes.metric)]
    } else {
      predicted.values <- predict(get_model(obj.axes, Env[[-(i - 3)]]),
                                  obj.axes@data)
      obj@variable.importance[(i - 3)] <- cor(predicted.values, o.predicted.values)
    }
  }

  # Variable importance normalization (%)
  if (sum(obj@variable.importance) == 0) {
    all.null <- TRUE
    for (i in seq_len(length(obj@variable.importance[1, ]))) {
      if (obj@variable.importance[1, i] != 0) {
        all.null <- FALSE
      }
    }
    if (all.null) {
      obj@variable.importance[1, ] <- 100/length(obj@variable.importance)
    } else {
      obj@variable.importance <- obj@variable.importance * 100
    }
  } else {
    obj@variable.importance <- obj@variable.importance/sum(obj@variable.importance) *
      100
  }
  row.names(obj@variable.importance) <- "Axes.evaluation"
  return(obj)
})
