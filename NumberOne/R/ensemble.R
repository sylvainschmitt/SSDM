#' @include Algorithm.Niche.Model.R
#' @importFrom raster raster stack reclassify
NULL

# Sum method used for ensemble function,
# Not documented because only an intern function
setMethod('sum', 'Algorithm.Niche.Model', function(x, ..., name = NULL, ensemble.metric = c('AUC'),
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
                variables.importance = x@variables.importance)
  smodel@data = smodel@data[-1,]
  smodel@variables.importance[1,] = 0

  # Name
  if (!is.null(name)) {name = paste0(name,'.')}
  smodel@name = paste0(name,class(x)[[1]],'.ensemble')

  # Datas, Projections, and Variables importance fusion
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
      weight.value = c(weight.value, ensemble.thresh[j])
    }
    weight.value = mean(weight.value)
    if (test) {
      if (weight) {
        smodel@projection = smodel@projection + models[[i]]@projection * weight.value
        smodel@variables.importance = smodel@variables.importance + models[[i]]@variables.importance * weight.value
        sweight = sweight + weight.value
      } else {
        smodel@projection = smodel@projection + models[[i]]@projection
        smodel@variables.importance = smodel@variables.importance + models[[i]]@variables.importance
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

      # Variables importance
      if (!is.numeric(sum(smodel@variables.importance))) {
        cat('Error variables importance is not numeric : \n')
        print(smodel@variables.importance)
        smodel@variables.importance = x@variables.importance
        smoadel@variable.importance[1,] = 0
      } else {
        if(is.nan(sum(smodel@variables.importance))) {
          cat('Error variables importance is NaN')
          smodel@variables.importance[1,] = (100/length(smodel@variables.importance))
        } else {
          if (sum(smodel@variables.importance) == 0) {
            all.null = T
            for(i in 1:length(smodel@variables.importance)) {if(smodel@variables.importance[1,i] != 0) {all.null = F}}
            if(all.null) {smodel@variables.importance[1,] = (100/length(smodel@variables.importance))} else {smodel@variables.importance = smodel@variables.importance * 100}
          }
          else {smodel@variables.importance = smodel@variables.importance / sum(smodel@variables.importance) * 100}
        }
      }

      # Evaluation
      smodel = evaluate(smodel, thresh = 1001)

      # Parameters
      smodel@parameters = x@parameters
      smodel@parameters$kept.model = kept.model

      return(smodel)
    }
})


#'Method to assemble different algorithms models in one ensemble model
#'
#'This is a method to assemble several algorithms models in one ensemble model.
#'It takes in inputs several S4 \linkS4class{Algorithm.Niche.Model} class
#'objects obtained with the \code{\link{Modelling}} function. It returns an S4
#'\linkS4class{Ensemble.Niche.Model} class object containing the habitat
#'suitability map, the binary map, and the uncertainty map based on the habitat
#'suitability map variance inter algorithms and all evaluation tables comming
#'with (model evaluation, algorithms evaluation, algorithms correlation matrix
#'and variables importance).
#'
#'@param x,... Algorithm.Niche.Model. Algortihms models to assemble.
#'@param name character. Optionnal name given to the final Ensemble.Niche.Model
#'  producted.
#'@param ensemble.metric character. Metric used to compute the selection among
#'  algorithms different models (see details below)
#'@param ensemble.thresh numeric. Threshold associated with the metric used to
#'  compute the selection.
#'@param weight logical. Choose if the model are weighted or not by the
#'  selection associated metrics mean.
#'@param thresh numeric. Binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'@param uncertainity logical. If false uncertainity mapping and algorithms
#'  correlation matrix is not computed.
#'
#'@details Ensemble metric (Metric used to compute the selection among
#'  algorithms different models) can be choosed among those : \describe{
#'  \item{AUC}{Area under the receiving operative curve (ROC)}
#'  \item{Kappa}{Kappa metric issued from the confusion matrix}
#'  \item{sensitivity}{Sensitivity metric issued from the confusion matrix}
#'  \item{specificity}{Specificity metric issued from the confusion matrix}
#'  \item{prop.correct}{Correct predicted occurences proportion issued from the
#'  confusion matrix} }
#'
#'@return an S4 \linkS4class{Ensemble.Niche.Model} Class object viewable with
#'  \code{\link{plot.model}} method
#'
#' @examples
#'\dontrun{
#' ensemble(GLM1, GLM2, GAM1, GAM2)
#'}
#'
#'@seealso \code{\link{Ensemble.Modelling}} for stack ensemble modelling with
#'  multiple algorithms
#'
#'@export
setMethod('ensemble', 'Algorithm.Niche.Model', function(x, ..., name = NULL,
                                                        ensemble.metric = c('AUC'), ensemble.thresh = c(0.75),
                                                        weight = T, thresh = 1001, uncertainity = T) {
  models = list(x, ...)
  enm = Ensemble.Niche.Model()

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
    cat('Projection, and variables importance computing...')
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
      enm@name = paste0(name,'Ensemble.Niche.Model')

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
      enm = evaluate(enm)
      cat('   done \n')

      # Axes evaluation
      cat('Axes evaluation...')
      enm@variables.importance = sum.algo.ensemble@variables.importance
      cat('   done \n')

      # Projections stack
      projections = stack()
      for (i in 1:length(algo.ensemble)) {
        projections = stack(projections, algo.ensemble[[i]]@projection)
        names(projections[[i]]) = algo.ensemble[[i]]@name
      }

      # Algorithms Correlation
      if (!(uncertainity)) {cat('Algorithm correlation computing is unactivated \n')}
      if (uncertainity && length(projections@layers) > 1) {
        cat('Algorithms correlation...')
        enm@algorithm.correlation = as.data.frame(layerStats(projections, 'pearson', na.rm = T)$`pearson correlation coefficient`)
        cat('   done \n')
      }

      # Uncertainity map
      if (!(uncertainity)) {cat('Uncertainty mapping is unactivated \n')}
      if (uncertainity && length(projections@layers) > 1) {
        cat('Uncertainity mapping...')
        enm@uncertainity = calc(projections, var)
        names(enm@uncertainity) = 'Uncertainity map'
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
