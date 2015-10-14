#### Modelling function #### ----
Modelling = function(algorithm,
                     # Modelling data input
                     Occurences, Env,
                     # Occurences reading
                     Xcol, Ycol, Pcol = NULL,
                     # Model creation
                     name = NULL,
                     # Pseudo-absences definition
                     PA = NULL, train.frac = 0.7,
                     # Evaluation parameters
                     thresh = 1001,
                     # Modelling parameters
                     ...) {
  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','RF')
  if(!(algorithm %in% available.algo)) {stop(algorithm,' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, RF')}
  
  # Empty Algorithm niche model object creation
  model = Algorithm.Niche.Model(algorithm)
  if (!is.null(name)) {name = paste0(name,'.')}
  model@name = paste0(name,algorithm,'.Niche.Model')
  
  cat('Data check ... \n')
  # Occurences data input test | Data frame needed
  if (is.matrix(Occurences)) {Occurences = data.frame(occurences)}
  if (!is.data.frame(Occurences)) {stop('Occurences data set is not a data frame or a matrix')}
  if ((Xcol %in% names(Occurences)) == F) {stop('X column is not well defined')}
  if ((Ycol %in% names(Occurences)) == F) {stop('Y column is not well defined')}
  if (is.null(Pcol)) {
    PO = T # Presence only
    cat('No presence column, presence-only data set is supposed.\n')
  } else if ((Pcol %in% names(Occurences)) == F) {stop('Presence column is not well defined')}
  if (!is.null(PA)) {PO = T}
  if (PO) {cat('Pseudo-absence selection will be computed.\n')}
  data = data.frame(X = Occurences[which(names(Occurences) == Xcol)], Y = Occurences[which(names(Occurences) == Ycol)])
  names(data) = c('X','Y')
  if (PO) {data$Presence = 1} else {data$Presence = Occurences[which(names(Occurences == Pcol))]}
  data$Train = F
  data$Train[sample.int(length(data$Presence), round(length(data$Presence)*train.frac))] = T
  
  # Environment data input test | RasterStack needed
  if (is.raster(Env)) {Env = stack(Env)}
  if (!inherits(Env, 'RasterStack')) {stop('Environment data set is not a raster or a raster stack')}
  cat('   done. \n\n')
  
  # Pseudo - absences selection
  cat('Pseudo absence selection... \n')
  model@data = data
  if (PO) {
    model = PA.select(model, Env, PA, train.frac)
    model@parameters['PA'] = T}
  model = data.values(model, Env)
  cat('   done. \n\n')
  
  # Projection
  cat('Model projection...')
  model = project(model, Env, ...)
  cat('   done. \n\n')
  
  # Evaluation
  cat('Model evaluation...\n')
  model = evaluate(model, thresh)
  cat('   done. \n\n')
  
  # Evaluation
  cat('Model axes contribution evaluation...\n')
  model = evaluate.axes(model, thresh, Env, ...)
  cat('   done. \n\n')
  
  return(model)
}

#### Ensemble Modelling function #### ----
Ensemble.Modelling = function(algorithms,
                              # Modelling data input
                              Occurences, Env,
                              # Occurences reading
                              Xcol, Ycol, Pcol = NULL,
                              # Model creation
                              rep = 1, name = NULL, save = F, directory = getwd(),
                              # Pseudo-absences definition
                              PA = NULL, train.frac = 0.7,
                              # Evaluation parameters
                              thresh = 1001, AUCthresh = 0.75, uncertainity = T,
                              # Modelling parameters
                              ...) {
  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','RF')
  for (i in 1:length(algorithms)) {
    if(!(algorithms[[i]] %in% available.algo)) {stop(algorithms[[i]],' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, RF')}}
  
  # Algorithms models creation
  cat('#### Algorithms models creation ##### \n\n')
  models = list()
  for (i in 1:length(algorithms)) {
    for (j in 1:rep) {
      model.name = paste0(algorithms[i],'.',j)
      cat('Modelling :', model.name, '\n\n')
      model = try(Modelling(algorithms[i], Occurences, Env, Xcol, Ycol, Pcol, name = NULL,
                            PA = NULL, train.frac = 0.7, thresh = 1001, ...))
      if (inherits(model, "try-error")) {cat(model)} else {models[model.name] = model}
      cat('\n\n')
    }
  }
  
  # Ensemble modelling
  cat('#### Ensemble modelling with algorithms models ##### \n\n')
  algo = list()
  for (i in 1:length(models)) {algo[[i]] = models[[i]]}
  algo['name'] = name
  algo['AUCthresh'] = AUCthresh
  algo['thresh'] = thresh
  algo['uncertainity'] = uncertainity
  enm = do.call(ensemble, algo)
  
  if(!is.null(enm)) {
    # Saving
    if(save) {
      cat('#### Saving ##### \n\n')
      save.enm(enm, directory)
    }
  }
  
  return(enm)
}
