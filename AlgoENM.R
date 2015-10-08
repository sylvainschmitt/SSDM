algorithm.modeling = function(data.algo, algo, name, thresh) {
  # modeling function for each algo
  algo.model = ENM()
  
  # GLM : Generalized Linear Model
  if(algo == 'GLM') {
    model = glm(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # GAM : Generalized Additive Model
  if(algo == 'GAM') {
    require(mgcv)
    model = gam(Presence ~ , data = data.algo)
    model = gam(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # MARS : Multivaraite Adaptative Regressions Splines
  if(algo == 'MARS') {
    require(earth)
    model = earth(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model)) # Projection issue with factor
    
  }
  
#   # MAXENT : Maximum Entropy
#   if(algo == 'GLM') {
#     require(dismo)
#     model = maxent(data.algo$Presence, data.algo)
#     algo.model@proj = stack(predict(Env, model))
#   }
  
  # ANN : Artificial Neural Network
  if(algo == 'ANN') {
    require(nnet)
    model = nnet(Presence ~ ., data = data.algo, size = 2)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # CTA : Classification Tree Analysis
  if(algo == 'CTA') {
    require(rpart)
    model = rpart(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # GBM : Generalized Boosted Model
  if(algo == 'GBM') {
    require(gbm)
    model = earth(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model)) # Projection issue with factor
  }
  
  # RF : Random Forest
  if(algo == 'RF') {
    require(randomForest)
    model = randomForest(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # Completing algo.model object
  algo.model@name = paste0(algo,'_',name)
  algo.model@data = data.algo
  threshold = optim.thresh(data.algo$Presence, model$fitted.values, thresh)
  threshold = mean(threshold$`max.sensitivity+specificity`)
  algo.model@eval = accuracy(data.algo$Presence, model$fitted.values, threshold)
  algo.model@axes.contrib = axes.contribution(data.algo, algo, thresh, algo.model@eval$AUC)
  
  return(model)
}

axes.contribution = function(data.algo, algo, name, thresh, AUC) {
  # Evaluate axes for each algo
  axes.contrib = data.frame(axes = character(), contribution = numeric())
  
  
  for (i in 2:length(data.algo)) {
    axes = names(data.algo)[i]
    
    # MARS : Multivaraite Adaptative Regressions Splines
    if(algo == 'MARS') {
      require(earth)
      model = earth(Presence ~ ., data = data.algo[-i])
      pred = fitted(model)
    }
    
    # RF : Random Forest
    if(algo == 'RF') {
      require(randomForest)
      model = randomForest(Presence ~ ., data = data.algo[-i])
      pred = predict(model, data.algo[-i])
    }
    
    # Computing AUC without each axes
    threshold = optim.thresh(data.algo$Presence, pred, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    value = AUC - accuracy(data.algo$Presence, pred, threshold)$AUC
    axes.contrib = rbind(axes.contrib, data.frame(axes = axes, contribution = value))
  }
  
  return(axes.contrib)
}


algorithm.modeling.loop = function(data, ...) {
  models = list()
  for (i in 1:length(data[grep('Run', names(data))])) {
    Run = data[grep(paste0('Run',i), names(data))] # selecting the run
    data.algo = data[which(Run == T),] # taking data for the algorithms
    data.algo = data.algo[-c(1:2)] # removing coordinates
    data.algo = data.algo[-grep('Run', names(data.algo))]
    model = algorithm.modeling(data.algo, ...)
    models[[i]] = model
  }
  return(models)
}