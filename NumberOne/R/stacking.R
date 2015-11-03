#' @include Ensemble.Niche.Model.R
#' @importFrom raster raster stack reclassify
NULL

#'Stack different ensemble models in one stack species model
#'
#'This is a function to stack several ensemble models in one stack species
#'model. It takes in inputs several S4 \linkS4class{Ensemble.Niche.Model} class
#'objects obtained with \code{\link{Ensemble.Modelling}} or
#'\code{\link{ensemble}} functions. It returns an S4
#'\linkS4class{Stack.Species.Ensemble.Niche.Model} class object containing the
#'local species richness map, and the uncertainty map based on the habitat
#'suitability map variance inter algorithms, all evaluation tables comming with
#'(model evaluation, algorithms evaluation, algorithms correlation matrix and
#'variables importance), and all associated ensemble models for each species
#'(see Ensemble.Modelling).
#'
#'@param enm,... character. Choice of the algorithm for the modelling (see
#'  details below).
#'@param name character. Optionnal name given to the final Ensemble.Niche.Model
#'  producted.
#'@param thresh numeric. binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'@param metric character. Method used to compute the binary map threshold (see
#'  details below.)
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'@param rep.B integer. If the method used to create the local species richness
#'  is random bernoulli (\strong{B}), it defines the number of repetition used
#'  to create random bernoulli binary maps for each species.
#'
#'@return an S4 \linkS4class{Stack.Species.Ensemble.Niche.Model} Class object
#'  viewable with \code{\link{plot}} method
#'
#'@details \strong{Metric :} choice of the metric used to compute binary map
#'  threshold and confusion matrix : \describe{ \item{"Kappa"}{maximizes the
#'  model Kappa value} \item{"TSS"}{\strong{True Skill Statistic} maximizes the
#'  sensitivity and specificity sum} \item{"CCR"}{maximizes the correct
#'  predicted observations proportion} \item{"SES"}{using the sensitivty
#'  specificity equality} \item{"LW"}{using the lowest occurence prediction
#'  probability} \item{"ROC"}{minimizing the distance between the ROC plot
#'  (receiving operative curve) and the upper left coin (1,1)} }
#' @examples
#'\dontrun{
#' stacking(Specie1.enm, Specie2.enm)
#'}
#'
#'@seealso \code{\link{Stack.Modelling}} for stack species ensemble modelling
#'  with multiple algorithms and multiples species
#'
#'@export
setMethod('stacking', 'Ensemble.Niche.Model', function(enm, ..., name = NULL, method = 'P',
                                                       metric = 'SES', thresh = 1001, rep.B = 1000) {
  enms = list(enm, ...)
  if (length(enms) < 2) {stop('You neeed more than one enm to do stackings')}
  names = c()
  for (i in 1:length(enms)) {if(enms[[i]]@name %in% names) {stop('Ensemble models can\'t have the same name, you need to rename one of ',enms[[i]]@name)} else {names = c(names, enms[[i]]@name)}}
  cat('Stack creation... \n')
  stack = Stack.Species.Ensemble.Niche.Model(diversity.map = reclassify(enm@projection[[1]], c(-Inf,Inf,0)),
                                             uncertainity = reclassify(enm@uncertainity, c(-Inf,Inf,NA)),
                                             parameters = enm@parameters)

  # Name
  cat('   naming...')
  if (!is.null(name)) {name = paste0(name,'.')}
  stack@name = paste0(name,'Stack.Species.Ensemble.Niche.Model')
  cat(' done. \n')

  # Diversity map
  cat('   diversity mapping...')
  # Useless datacheck to prevent bugs to remove after debugging
  for (i in 1:length(enms)) {
    if(!inherits(enms[[i]]@projection, 'RasterLayer')){
      cat('Error', enms[[i]]@name, 'is not a raster but a', class(enms[[i]]@projection)[1], '.\nIt will be removed for the stacking')
      enms[[i]] = NULL
    }
  }
  if (method == 'P') {
    cat('\n Local species richness coomputed by summing individual probabilities. \n')
    for (i in 1:length(enms)) {stack@diversity.map = stack@diversity.map + enms[[i]]@projection}
  }
  if (method == 'T') {
    cat('\n Local species richness coomputed by thresholding and then summing. \n')
    for (i in 1:length(enms)) {
      enms[[i]] = evaluate(enms[[i]], thresh = thresh, metric = metric)
      stack@diversity.map = stack@diversity.map +
        reclassify(enms[[i]]@projection,
                   c(-Inf,enms[[1]]@evaluation$threshold,0, enms[[i]]@evaluation$threshold,Inf,1))}
  }
  if (method == 'B') {
    cat('\n Local species richness coomputed by drawing repeatedly from a Bernoulli distribution. \n')
    # rbinom(lengths(enms), 1000 trials, enms.proba)
    proba = stack()
    for (i in 1:length(enms)) {proba = stack(proba, enms[[i]]@projection)}
    diversity.map = calc(proba, fun = function(...) {
      x = c(...)
      x[is.na(x)] = 0
      return(rbinom(lengths(x), rep.B, x))},
      forcefun = T)
    stack@diversity.map = sum(diversity.map) / length(enms) / rep.B
  }
  names(stack@diversity.map) = 'diversity'
  cat(' done. \n')

  # Uncertainity map
  cat('   uncertainity mapping...')
  uncertainities = stack()
  for (i in 1:length(enms)) {
    a = try(enms[[i]]@uncertainity)
    if (inherits(a, 'try-error')) {cat('Ensemble model',enms[[i]]@name,'uncertinity map not computed')} else {uncertainities = stack(uncertainities, a)}
  }
  a = try(calc(uncertainities, mean))
  if (inherits(a, 'try-error')) {cat('No uncertainity map to do uncertainity mapping')
  } else {
    stack@uncertainity = a
    names(stack@uncertainity) = 'uncertainity'
  }

  cat(' done. \n')

  # Evaluation
  cat('   evaluating...')
  stack@evaluation = enm@evaluation
  for (i in 2:length(enms)) {stack@evaluation = rbind(stack@evaluation, enms[[i]]@evaluation)}
  a = stack@evaluation[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in 1:length(stack@evaluation)) {a[i] = c(mean(stack@evaluation[,i], na.rm = T), sd(stack@evaluation[,i], na.rm = T))}
  stack@evaluation = a
  cat(' done. \n')

  # Variables Importance
  cat('   comparing variables importnace...')
  stack@variables.importance = enm@variables.importance
  for (i in 2:length(enms)) {
    a  = try(rbind(stack@variables.importance, enms[[i]]@variables.importance))
    if (inherits(a, 'try-error')) {cat(a)} else {stack@variables.importance = a}
  }
  a = stack@variables.importance[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in 1:length(stack@variables.importance)) {a[i] = c(mean(stack@variables.importance[,i]), sd(stack@variables.importance[,i]))}
  stack@variables.importance = a
  cat(' done. \n')

  # Algorithm Correlation
  cat('   comparing algorithms correlation...')

  algo = c() # Listing all algorithms presents in enms and renaming enms row and columns
  for (i in 1:length(enms)) {
    if(length(enms[[i]]@algorithm.correlation) == 0) {cat('\n', enms[[i]]@name,'algorithms correlation has not been computed. \n')} else {
      for (j in 1:length(enms[[i]]@algorithm.correlation)) {
        if (length(strsplit(names(enms[[i]]@algorithm.correlation)[j], '.', fixed = T)[[1]]) > 1){
          names(enms[[i]]@algorithm.correlation)[j] = strsplit(names(enms[[i]]@algorithm.correlation)[j], '.', fixed = T)[[1]][2]
          row.names(enms[[i]]@algorithm.correlation)[j] = strsplit(row.names(enms[[i]]@algorithm.correlation)[j], '.', fixed = T)[[1]][2]
        }
        if(!(names(enms[[i]]@algorithm.correlation)[j] %in% algo)) {algo = c(algo, names(enms[[i]]@algorithm.correlation)[j])}
      }
    }
  }
  mcorr = data.frame(matrix(nrow = length(algo), ncol = length(algo)))
  names(mcorr) = algo
  row.names(mcorr) = algo
  if(length(algo) > 0) {
    for (i in 1:length(algo)) {
      for (j in 1:length(algo)) {
        if(i > j) {
          corr = c()
          for (k in 1:length(enms)) {
            if(length(enms[[k]]@algorithm.correlation) != 0) {
              row = which(row.names(enms[[k]]@algorithm.correlation) == row.names(mcorr)[j])
              col = which(names(enms[[k]]@algorithm.correlation) == names(mcorr)[i])
              if(length(row) > 0 && length(col) > 0) {corr = c(corr, enms[[k]]@algorithm.correlation[row,col])}
            }
            mcorr[i,j] = mean(corr, na.rm = T)
          }
        }
      }
    }
  }
  stack@algorithm.correlation = mcorr
  cat(' done. \n')

  # Algorithm Evaluation
  cat('   comparing algorithms evaluation')
  stack@algorithm.evaluation = enm@algorithm.evaluation
  for (i in 2:length(enms)) {stack@algorithm.evaluation = rbind(stack@algorithm.evaluation, enms[[i]]@algorithm.evaluation)}
  stack@algorithm.evaluation$algo = 'algo'
  for (i in 1:length(row.names(stack@algorithm.evaluation))) {stack@algorithm.evaluation$algo[i] = strsplit(row.names(stack@algorithm.evaluation),'.', fixed = T)[[i]][2]}
  stack@algorithm.evaluation = aggregate.data.frame(stack@algorithm.evaluation, by = list(stack@algorithm.evaluation[,which(names(stack@algorithm.evaluation) == 'algo')]), FUN = mean)
  row.names(stack@algorithm.evaluation) = stack@algorithm.evaluation$Group.1
  stack@algorithm.evaluation = stack@algorithm.evaluation[-1]
  stack@algorithm.evaluation = stack@algorithm.evaluation[-which(names(stack@algorithm.evaluation) == 'algo')]
  cat(' done. \n')

  # ENMS
  for (i in 1:length(enms)) {stack@enms[enms[[i]]@name] = enms[[i]]}

  # Parameters
  stack@parameters$method = method
  if (method == 'B') {stack@parameters$rep.B = rep.B}

  return(stack)})
