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