predict.algo = function(model, x){
  x= as.data.frame(x)
  for (i in 1:length(Env@layers)) {
    if(Env[[i]]@data@isfactor) {
      x[,i] = as.factor(Env[[i]]@data@attributes[[1]]$ID[x[,i]])
      levels(x[,i]) = Env[[i]]@data@attributes[[1]]$ID
    }
  }
  return(predict(model, x))
}
make.gam.formula = function(data.algo) {
  formula = "Presence ~"
  for (i in 2:length(data.algo)) {
    var = names(data.algo[i])
    if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
    if (!is.factor(data.algo[,i])) {formula = paste0(formula,' + s(',var,')')}
  }
  return(formula)
}
axes.contribution = function(data.algo, algo, name, thresh, AUC) {
  # Evaluate axes for each algo
  axes.contrib = data.frame(axes = character(), contribution = numeric())
  for (i in 2:length(data.algo)) {
    axes = names(data.algo)[i]
    
    
    # Computing AUC without each axes
    threshold = optim.thresh(data.algo$Presence, pred, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    value = AUC - accuracy(data.algo$Presence, pred, threshold)$AUC
    axes.contrib = rbind(axes.contrib, data.frame(axes = axes, contribution = value))
  }
  
  return(axes.contrib)
}

algorithm.modeling = function(data.algo, Env, # data
                              algo, name = '', # choices
                              test = 'AIC', epsilon = 1e-08, maxit = 500, # GLM options
                              final.leave = 1, cv = 3, # CTA options
                              thresh.shrink = 1e-03, trees = 2500, # GBM options
                              thresh # threshold computing precision
                              ) {
  # Modeling function for each algo
  algo.model = ENM()
  
  # GLM : Generalized Linear Model
  if(algo == 'GLM') {
    model = glm(Presence ~ ., data = data.algo, 
                test = test, control = glm.control(epsilon = epsilon, maxit = maxit))
    algo.model@proj = stack(predict(Env, model))
  }
  
  # GAM : Generalized Additive Model
  if(algo == 'GAM') {
    require(mgcv)
    formula = make.gam.formula(data.algo)
    model = gam(formula(formula), data = data.algo)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # MARS : Multivaraite Adaptative Regressions Splines
  if(algo == 'MARS') {
    require(earth)
    model = earth(Presence ~ ., data = data.algo, degree = 2)
    algo.model@proj = stack(predict(Env, model, fun = predict.algo))
  }
  
  # MAXENT : Maximum Entropy
  if(algo == 'MAXENT') {
    require(dismo)
    data.algo.presence = data.algo[-which(data.algo$Presence == 0),]
    p = data.algo.presence[c(1:2)]

    factors = c()
    x = data.algo.presence[-c(1:3)]
    for (i in 1:length(names(x))) {
      if(is.factor(x[,i])) {factors = c(factors, names(x)[i])}
    }
    model = maxent(x = Env, p = p, factors = factors)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # ANN : Artificial Neural Network
  if(algo == 'ANN') {
    require(nnet)
    model = nnet(Presence ~ ., data = data.algo, size = 2, maxit = maxit)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # CTA : Classification Tree Analysis
  if(algo == 'CTA') {
    require(rpart)
    model = rpart(Presence ~ ., data = data.algo, method = 'class', control = rpart.control(minbucket = final.leave, xval = cv))
    algo.model@proj = stack(predict(Env, model))
  }
  
  # GBM : Generalized Boosted Model
  if(algo == 'GBM') {
    require(gbm)
    model = gbm(Presence ~ ., data = data.algo,
                distribution = 'bernoulli', n.minobsinnode = final.leave,
                shrinkage = thresh.shrink, bag.fraction = 0.5, 
                train.fraction = 1, cv.folds = cv, n.trees = trees)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # RF : Random Forest
  if(algo == 'RF') {
    require(randomForest)
    model = randomForest(Presence ~ ., data = data.algo, do.classif = TRUE, ntree = trees, nodesize = final.leave, maxnodes = NULL)
    algo.model@proj = stack(predict(Env, model))
  }
  
  # SVM : Support Vector Machine
  if(algo == 'SVM') {
    require(e1071)
    model = svm(Presence ~ ., data = data.algo, type = 'eps-regression', gamma = 1/(length(data.algo)-1), kernel = 'radial', epsilon = epsilon, cross = cv)
    algo.model@proj = stack(predict(Env, model, fun = predict.algo))
  }
  
  # FDA : Flexible Discriminant Analysis
  if(algo == 'FDA') {
    require(mda)
    model = fda(Presence ~ ., data = data.algo)
    algo.model@proj = stack(predict(Env, model, fun = predict.algo))
  }
  
  # Rescaling projection
  proj = algo.model@proj[[1]]
  if(proj@data@min < 0) {proj = reclassify(proj, c(-Inf,0,0))}
  proj = proj / proj@data@max
  algo.model@proj[[1]] = proj
  
  # Completing algo.model object
#   algo.model@name = paste0(algo,'_',name)
#   algo.model@data = data.algo
#   require(SDMTools)
#   threshold = optim.thresh(data.algo$Presence, model$fitted.values, thresh)
#   threshold = mean(threshold$`max.sensitivity+specificity`)
#   algo.model@eval = accuracy(data.algo$Presence, model$fitted.values, threshold)
#   algo.model@axes.contrib = axes.contribution(data.algo, algo, thresh, algo.model@eval$AUC)
  
  return(algo.model)
}

algorithm.modeling.loop = function(data, ...) {
  models = list()
  for (i in 1:length(data[grep('Run', names(data))])) {
    Run = data[grep(paste0('Run',i), names(data))] # selecting the run
    data.algo = data[which(Run == T),] # taking data for the algorithms
    data.algo = data.algo[-grep('Run', names(data.algo))]
    data.algo = data.algo[-c(1:2)] # removing coordinates
    model = algorithm.modeling(data.algo, ...)
    models[[i]] = model
  }
  return(models)
}