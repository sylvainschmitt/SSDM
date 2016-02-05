#' @include Algorithm.SDM.R checkargs.R
#' @importFrom raster raster stack reclassify
NULL

#'Methods to assemble multiple algorithms in an ensemble SDM
#'
#'This is a method to assemble several algorithms in an ensemble SDM. The
#'function takes as inputs several S4 \linkS4class{Algorithm.SDM} class objects
#'obtained with the \code{\link{modelling}} function. The function returns an S4
#'\linkS4class{Ensemble.SDM} class object containing the habitat suitability
#'map, the binary map, and the uncertainty map (based on the between-algorithm
#'variance) and the associated evaluation tables (model evaluation,
#'algorithm evaluation, algorithm correlation matrix and variable importance).
#'
#'@param x,... SDMs. SDMs to be assembled.
#'@param name character. Optional name given to the final Ensemble.SDM produced
#'  (by default 'Ensemble.SDM').
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'@param weight logical. Choose whether or not you want the SDMs to be weighted
#'  using the selection metric or, alternatively, the mean of the selection
#'  metrics.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param uncertainty logical. If set to true, generates an uncertainty map and an algorithm
#'  correlation matrix.
#'@param verbose logical. If set to true, allows the function to print text in the
#'  console.
#'@param GUI,format,na.rm  logical. Don't take those arguments into account
#'  (parameters for the user interface and sum function).
#'
#'@details ensemble.metric (metric(s) used to select the best SDMs that will be
#'  included in the ensemble SDM) can be chosen from among: \describe{
#'  \item{AUC}{Area under the receiving operating characteristic (ROC) curve}
#'  \item{Kappa}{Kappa from the confusion matrix} \item{sensitivity}{Sensitivity
#'  from the confusion matrix} \item{specificity}{Specificity from the confusion
#'  matrix} \item{prop.correct}{Proportion of correctly predicted occurrences
#'  from the confusion matrix} }
#'
#'@return an S4 \linkS4class{Ensemble.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occurrences = subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#'
#' # ensemble SDM building
#' CTA = modelling('CTA', Occurrences, Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE')
#' SVM = modelling('SVM', Occurrences, Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE')
#' ESDM = ensemble(CTA, SVM, ensemble.thresh = c(0.6))
#'
#' # Results plotting
#' plot(ESDM)
#' }
#'
#'@seealso \code{\link{ensemble_modelling}} to build an ensemble SDM from
#'  multiple algorithms.
#'
#'@name ensemble
#'
#'@export
setGeneric('ensemble',
           function(x, ..., name = NULL,ensemble.metric = c('AUC'), ensemble.thresh = c(0.75),
                    weight = T, thresh = 1001, uncertainty = T, verbose = T, GUI = F) {return(standardGeneric('ensemble'))})

#' @rdname ensemble
#' @export
setMethod('ensemble', 'Algorithm.SDM',
          function(x, ..., name = NULL, ensemble.metric = c('AUC'), ensemble.thresh = c(0.75),
                   weight = T, thresh = 1001, uncertainty = T, verbose = T, GUI = F) {
  # Check agruments
  .checkargs(name = name, ensemble.metric =ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, thresh = thresh, uncertainty = uncertainty, verbose = verbose, GUI = GUI)

  models = list(x, ...)
  enm = Ensemble.SDM()

  # Algorithm ensemble model creation
  cat('Creation of one ensemble niche model by algorithm...')
  algo.ensemble = list()
  while(length(models) > 0) {
    type.model = list()
    type = class(models[[1]])[[1]]
    rm = {}
    for (i in 1:length(models)) {
      if (inherits(models[[i]], type)) {
        type.model[(length(type.model)+1)] = models[[i]]
        rm = c(rm, i)
      }
    }
    if (length(rm) > 0) {
      for (i in 1:length(rm)) {
        models[[rm[i]]] = NULL
        rm = rm - 1
      }
    }
    type.model['name'] = name
    type.model[['ensemble.metric']] = ensemble.metric
    type.model[['ensemble.thresh']] = ensemble.thresh
    type.model['weight'] = weight
    type.model['thresh'] = thresh
    type.model['format'] = T
    type.model['verbose'] = F
    algo.ensemble[type] = do.call(sum, type.model)
  }
  cat('   done. \n')

  if (length(algo.ensemble) < 1) {
    cat('No model were kept with this threshold, Null is returned. \n')
    return(NULL)
  } else {

    # Sum of algorithm ensemble
    cat('Projection, evaluation and variable importance computing...')
    algo.list = list()
    for (i in 1:length(algo.ensemble)) {algo.list[[i]] = algo.ensemble[[i]]}
    algo.list['name'] = 'sum'
    algo.list[['ensemble.metric']] = ensemble.metric
    algo.list[['ensemble.thresh']] = ensemble.thresh
    algo.list['weight'] = weight
    algo.list['thresh'] = thresh
    algo.list['format'] = F
    sum.algo.ensemble = do.call(sum, algo.list)
    if (length(sum.algo.ensemble) < 1) {
      return(NULL)
    } else {

      # Name
      if (!is.null(name)) {name = paste0(name,'.')} else {name = 'Specie.'}
      enm@name = paste0(name,'Ensemble.SDM')

      # Projection
      enm@projection = sum.algo.ensemble@projection
      cat('   done \n')

      # Data
      enm@data = algo.ensemble[[1]]@data
      if (length(algo.ensemble) > 1) {
        for (i in 2:length(algo.ensemble)) {
          enm@data = rbind(enm@data, algo.ensemble[[i]]@data)
        }
      }

      # Evaluation
      cat('Model evaluation...')
      enm@evaluation = sum.algo.ensemble@evaluation
      cat('   done \n')

      # Axes evaluation
      cat('Axes evaluation...')
      enm@variable.importance = sum.algo.ensemble@variable.importance
      cat('   done \n')

      # Projections stack
      projections = stack()
      for (i in 1:length(algo.ensemble)) {
        projections = stack(projections, algo.ensemble[[i]]@projection)
        names(projections[[i]]) = algo.ensemble[[i]]@name
      }

      # Algorithms Correlation
      if (!(uncertainty)) {cat('Algorithm correlation computing is unactivated \n')}
      if (uncertainty && length(projections@layers) > 1) {
        cat('Algorithms correlation...')
        enm@algorithm.correlation = as.data.frame(layerStats(projections, 'pearson', na.rm = T)$`pearson correlation coefficient`)
        cat('   done \n')
      }

      # uncertainty map
      if (!(uncertainty)) {cat('Uncertainty mapping is unactivated \n')}
      if (uncertainty && length(projections@layers) > 1) {
        cat('uncertainty mapping...')
        enm@uncertainty = calc(projections, var)
        names(enm@uncertainty) = 'uncertainty map'
        cat('   done \n')
      }

      # Algorithms Evaluation
      cat('Algorithms evaluation...')
      enm@algorithm.evaluation = algo.ensemble[[1]]@evaluation
      row.names(enm@algorithm.evaluation)[1] = algo.ensemble[[1]]@name
      if (length(algo.ensemble) > 1) {
        for (i in 2:length(algo.ensemble)) {
          enm@algorithm.evaluation = rbind(enm@algorithm.evaluation, algo.ensemble[[i]]@evaluation)
          row.names(enm@algorithm.evaluation)[i] = algo.ensemble[[i]]@name
        }
      }
      enm@algorithm.evaluation$kept.model = algo.ensemble[[1]]@parameters$kept.model
      if (length(algo.ensemble) > 1) {
        for (i in 2:length(algo.ensemble)) {
          enm@algorithm.evaluation$kept.model[[i]] = algo.ensemble[[i]]@parameters$kept.model
        }

      }

      # Parameters
      enm@parameters = algo.ensemble[[1]]@parameters
      text.ensemble.metric = character()
      text.ensemble.thresh = character()
      for (i in 1:length(ensemble.metric)) {
        text.ensemble.metric = paste0(text.ensemble.metric,'.',ensemble.metric[i])
        text.ensemble.thresh = paste0(text.ensemble.thresh,'|',ensemble.thresh[i])
      }
      enm@parameters$ensemble.metric = text.ensemble.metric
      enm@parameters$ensemble.thresh = text.ensemble.thresh
      enm@parameters$weight = weight
    }

    cat('   done \n')

    return(enm)}})

#' @rdname ensemble
#' @export
setMethod('sum', 'Algorithm.SDM', function(x, ..., name = NULL, ensemble.metric = c('AUC'),
                                                   ensemble.thresh = c(0.75), weight = T,
                                                   thresh = 1001, format = T, verbose = T, na.rm = F) {
  models = list(x, ...)
  if (length(ensemble.metric) != length(ensemble.thresh)) {stop('You must have the same number of metrics and associated thresholds in models assembling step (see ensemble.metric and ensemble.thresh)')}
  if(format) {
    for(i in 1:length(models)) {
      if(!inherits(models[[i]], class(x)[[1]])) {
        stop('You can only sum models from the same algorithm')
      }
    }
  }

  smodel =  new(class(x)[[1]],
                projection = reclassify(x@projection[[1]], c(-Inf,Inf,0)),
                data = x@data[1,],
                variable.importance = x@variable.importance,
                evaluation = x@evaluation)
  smodel@data = smodel@data[-1,]
  smodel@variable.importance[1,] = 0
  smodel@evaluation[1,] = 0

  # Name
  if (!is.null(name)) {name = paste0(name,'.')}
  smodel@name = paste0(name,class(x)[[1]],'.ensemble')

  # Datas, Projections, Evaluation and variable importance fusion
  sweight = 0
  kept.model = 0
  for (i in 1:length(models)) {
    # Assembling selection test depending on parameters
    test = T
    weight.value = c()
    for (j in 1:length(ensemble.metric)) {
      if(models[[i]]@evaluation[,which(names(models[[i]]@evaluation) == ensemble.metric[j])] < ensemble.thresh[j]) {
        test = F
      }
      weight.value = c(weight.value, models[[i]]@evaluation[,which(names(models[[i]]@evaluation) == ensemble.metric[j])])
    }
    weight.value = mean(weight.value)
    if (test) {
      if (weight) {
        smodel@projection = smodel@projection + models[[i]]@projection * weight.value
        smodel@variable.importance = smodel@variable.importance + models[[i]]@variable.importance * weight.value
        smodel@evaluation = smodel@evaluation + models[[i]]@evaluation * weight.value
        sweight = sweight + weight.value
      } else {
        smodel@projection = smodel@projection + models[[i]]@projection
        smodel@variable.importance = smodel@variable.importance + models[[i]]@variable.importance
        smodel@evaluation = smodel@evaluation + models[[i]]@evaluation
        sweight = sweight + 1
      }
      smodel@data = rbind(smodel@data, models[[i]]@data)
      kept.model = kept.model + 1
    }
  }

  # Return NULL if any model is kept
  if( kept.model == 0) {
    if (verbose) {cat('No model were kept with this threshold, Null is return. \n')}
    return(NULL)} else {

      smodel@projection = smodel@projection / sweight
      names(smodel@projection) = 'Probability'
      smodel@evaluation = smodel@evaluation / sweight

      # variable importance
      if (!is.numeric(sum(smodel@variable.importance))) {
        cat('Error variable importance is not numeric : \n')
        print(smodel@variable.importance)
        smodel@variable.importance = x@variable.importance
        smoadel@variable.importance[1,] = 0
      } else {
        if(is.nan(sum(smodel@variable.importance))) {
          cat('Error variable importance is NaN')
          smodel@variable.importance[1,] = (100/length(smodel@variable.importance))
        } else {
          if (sum(smodel@variable.importance) == 0) {
            all.null = T
            for(i in 1:length(smodel@variable.importance)) {if(smodel@variable.importance[1,i] != 0) {all.null = F}}
            if(all.null) {smodel@variable.importance[1,] = (100/length(smodel@variable.importance))} else {smodel@variable.importance = smodel@variable.importance * 100}
          }
          else {smodel@variable.importance = smodel@variable.importance / sum(smodel@variable.importance) * 100}
        }
      }

      # Parameters
      smodel@parameters = x@parameters
      smodel@parameters$kept.model = kept.model

      return(smodel)
    }
})
