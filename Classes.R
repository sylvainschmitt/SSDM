##### Libraries ##### ----
library(raster)
library(rgdal)
library(methods)
library(SDMTools)
library(mgcv)
library(earth)
library(rpart)
library(gbm)
library(randomForest)
library(dismo)
library(nnet)
library(e1071)
library(shiny)
library(shinydashboard)
library(gplots)

##### New generics ##### ----
setGeneric('evaluate', function(obj, thresh = 1001) {return(standardGeneric('evaluate'))})
setGeneric('get_PA', function(obj) {return(standardGeneric('get_PA'))})
setGeneric('PA.select', function(obj, Env, ...) {return(standardGeneric('PA.select'))})
setGeneric('data.values', function(obj, Env, na.rm = T) {return(standardGeneric('data.values'))})
setGeneric('get_model', function(obj, ...) {return(standardGeneric('get_model'))})
setGeneric('project', function(obj, Env, ...) {return(standardGeneric('project'))})
setGeneric('evaluate.axes', function(obj, thresh = 1001, Env, ...) {return(standardGeneric('evaluate.axes'))})
setGeneric('ensemble', function(x, ..., name = NULL, AUCthresh = 0.75, thresh = 1001, uncertainity = T) {return(standardGeneric('ensemble'))})
setGeneric('save.enm', function (enm, directory = getwd()) {return(standardGeneric('save.enm'))})
setGeneric('load.enm', function (directory = getwd()) {return(standardGeneric('load.enm'))})
setGeneric('stacking', function(enm, ...) {return(standardGeneric('stacking'))})

##### Niche Model Class ##### -----

# 1 - Class definition #
setClass('Niche.Model',
         representation(name = 'character', 
                        projection = 'RasterStack', 
                        evaluation = 'data.frame', 
                        variables.importance = 'data.frame',
                        data = 'data.frame',
                        parameters = 'list'),
         prototype(name = character(), 
                   projection = stack(), 
                   evaluation = data.frame(), 
                   variables.importance = data.frame(),
                   data = data.frame(),
                   parameters = list()))

# 2 - Methods definition #
setMethod("evaluate", "Niche.Model", function(obj, thresh = 1001) {
  data = obj@data[which(!obj@data$Train),]
  predicted.values = extract(obj@projection[[1]], data[c('X','Y')])
  predicted.values[which(is.na(predicted.values))] = 0
  threshold = optim.thresh(data$Presence, predicted.values, thresh)
  threshold = mean(threshold$`max.sensitivity+specificity`)
  obj@evaluation = accuracy(data$Presence, predicted.values, threshold)
  row.names(obj@evaluation) = "Evaluation"
  return(obj)})

setMethod('print', 'Niche.Model', function(x, ...) {
  cat('Object of class :', class(x)[1],'\n')
  cat('Name :', x@name, '\n')
  cat('Projections : ',names(x@projection),'\n')
  print(x@evaluation)
  print(x@variables.importance)
  if(inherits(x, 'Ensemble.Niche.Model')) {
    cat('Uncertinity map :', names(x@uncertainity),'\n')
    print(x@algorithm.evaluation)
    print(x@algorithm.correlation)
  }
})

setMethod('plot', 'Niche.Model', function(x, y, ...) {
  if (inherits(x, 'Algorithm.Niche.Model')) {full = F} else {full = T}
  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(disable = T),
    dashboardBody(
      fluidRow(
        tabBox(title = 'Maps',
               tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                         plotOutput('probability', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)), title = 'Probability'),
               if(full) {tabPanel(plotOutput('niche'), title = 'Niche')},
               if(full) {tabPanel(plotOutput('uncertainity'), title = 'Uncertainity')}
        ),
        tabBox(title = 'Variables importance',
               tabPanel(plotOutput('varimp.barplot'), title = 'Barplot'), 
               tabPanel(tableOutput('varimp.table'), title = 'Table')
        )
      ),
      if(full) {
        fluidRow(
          tabBox(title = 'Model evaluation',
                 tabPanel(plotOutput('evaluation.barplot'), title = 'Barplot'), 
                 tabPanel(tableOutput('evaluation.table'), title = 'Table')
          ),
          if (length(x@algorithm.correlation) > 0) {
            tabBox(title = 'Algorithms correlation',
                   tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'), 
                   tabPanel(tableOutput('algo.corr.table'), title = 'Table')
            )}
        )}
    )
  )
  
  
  server <- function(input, output) {
    set.seed(122)
    
    if(full) {
      for (i in 1:length(row.names(x@algorithm.evaluation))) {row.names(x@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@algorithm.evaluation)[i]), '.Niche')[[1]][1]}
      if (length(x@algorithm.correlation) > 0) {
        for (i in 1:length(row.names(x@algorithm.correlation))) {row.names(x@algorithm.correlation)[i] = strsplit(as.character(row.names(x@algorithm.correlation)[i]), '.', fixed = T)[[1]][2]}
        x@algorithm.correlation[upper.tri(x@algorithm.correlation, diag = T)] = NA
        names(x@algorithm.correlation) = row.names(x@algorithm.correlation)
      }
    }
      # Maps
      
      # Single zoomable plot
      ranges <- reactiveValues(x = NULL, y = NULL)
      # When a double-click happens, check if there's a brush on the plot.
      # If so, zoom to the brush bounds; if not, reset the zoom.
      observeEvent(input$plot1_dblclick, {
        brush <- input$plot1_brush
        if (!is.null(brush)) {
          ranges$x <- c(brush$xmin, brush$xmax)
          ranges$y <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges$x <- NULL
          ranges$y <- NULL
        }
      })
      observeEvent(input$unzoom, {
        ranges$x <- NULL
        ranges$y <- NULL
      })
      
      output$probability <- renderPlot({
        if (!is.null(ranges$x)) {proba.map = crop(x@projection[[1]], c(ranges$x, ranges$y))} else {proba.map = x@projection[[1]]}
        plot(proba.map, 
             main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)),
             xlab = 'Longitude (°)',
             ylab = 'Latitude (°)')
        points(x@data$X[which(x@data$Presence == 1)], 
               x@data$Y[which(x@data$Presence == 1)], 
               pch = 16, cex = 0.7)
      })
      output$niche <- renderPlot({
        if (!is.null(ranges$x)) {niche.map = crop(x@projection[[2]], c(ranges$x, ranges$y))} else {niche.map = x@projection[[2]]}
        plot(niche.map, main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)))})
      if(full) {
        output$uncertainity <- renderPlot({
          if (!is.null(ranges$x)) {uncert.map = crop(x@uncertainity, c(ranges$x, ranges$y))} else {uncert.map = x@uncertainity}
          plot(uncert.map, main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)))})
        # Evaluation
        output$evaluation.barplot <- renderPlot({
          evaluation = x@algorithm.evaluation
          evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
          table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
          barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
          legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
        })
        output$evaluation.table <- renderTable({x@algorithm.evaluation[c(2,4:8)]})
        if (length(x@algorithm.correlation) > 0) {
          # Algorithms correlation
          output$algo.corr.table <- renderTable({x@algorithm.correlation})
          output$algo.corr.heatmap <- renderPlot({
            m <- as.matrix(x@algorithm.correlation)
            heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                      cellnote = round(m,3), notecol = "black", notecex = 2,
                      trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
          })
        }
      }
      # Variable importance
      output$varimp.barplot <- renderPlot({
        varimp = as.data.frame(t(x@variables.importance[-1]))
        names(varimp) = 'Axes.evaluation'
        barplot(varimp$Axes.evaluation, names.arg = row.names(varimp), las = 2)
      })
      output$varimp.table <- renderTable({x@variables.importance[-1]})
    }
    
    shinyApp(ui, server)
  })

##### Algorithm Niche Model Class ##### -----

# 1 - Class definition #
setClass('Algorithm.Niche.Model',
         contains = 'Niche.Model')

# 2 - Class creation function #
Algorithm.Niche.Model <- function(algorithm = 'Algorithm',
                                  name = character(), 
                                  projection = stack(), 
                                  evaluation = data.frame(), 
                                  variables.importance = data.frame(),
                                  data = data.frame(),
                                  parameters = list()) {
  object.class = paste0(algorithm,'.Niche.Model')
  return(new(object.class, name = name, projection = projection, evaluation = evaluation, variables.importance = variables.importance, data = data, parameters = list()))
}

# 3 - Methods definition #
setMethod('get_PA', "Algorithm.Niche.Model", function(obj) {return(obj)})

setMethod('PA.select', "Algorithm.Niche.Model", function(obj, Env, PA = NULL, train.frac = 0.7) {
  if (is.null(PA)) {PA = get_PA(obj)}
  
  # Mask defining
  if (PA$strat == '2nd') {
    cat('   second far selection \n')
    circles = list()
    for (i in 1:length(obj@data$X)) {
      x = obj@data$X[i]
      y = obj@data$Y[i]
      pts = seq(0, 2 * pi, length.out = 100)
      xy = cbind(x + 2/60 * sin(pts), y + 2/60 * cos(pts))
      circle = Polygon(xy)
      circles[i] = circle
    }
    sc= SpatialPolygons(list(Polygons(circles, 'Circles')))
    Mask = mask(Env[[1]], sc)
  } else {
    cat('   random selection \n')
    Mask = Env[[1]]}
  
  # Pseudo-Absences selection
  data.PA = data.frame(matrix(nrow = 0, ncol = 2))
  names(data.PA) = c('X','Y')
  if(PA$nb < 100) {nb = PA$nb*PA$nb} else {nb = 1000}
  while (length(data.PA[,1]) < PA$nb) {
    X = runif(nb, min = bbox(Mask)[1,1], max = bbox(Mask)[1,2])
    Y = runif(nb, min = bbox(Mask)[2,1],max = bbox(Mask)[2,2])
    points = data.frame(X = X, Y = Y)
    check = extract(Mask, points)
    points = points[-which(is.na(check)),]
    data.PA = rbind(data.PA, points)
  }
  data.PA = data.PA[1:PA$nb,]
  data.PA$Presence = 0
  data.PA$Train = F
  data.PA$Train[sample.int(length(data.PA$Presence), round(length(data.PA$Presence)*train.frac))] = T
  obj@data = rbind(obj@data, data.PA)
  
  return(obj)})

setMethod('data.values', "Algorithm.Niche.Model", function(obj, Env, na.rm = T) {
  values = data.frame(extract(Env, cbind(obj@data$X, obj@data$Y)))
  
  # Categorical variables as factor
  for (i in 1:length(Env@layers)) {
    if(Env[[i]]@data@isfactor) {
      col = which(names(values) == Env[[i]]@data@names)
      values[,col] = as.factor(values[,col])
      levels(values[,col]) = Env[[i]]@data@attributes[[1]]$ID
    }
  }
  
  # Tables binding
  obj@data = cbind(obj@data, values)
  
  # NAs removing
  if(na.rm) {
    for (i in 1:length(Env@layers)) {
      if(length(which(is.na(obj@data[i+3]))) > 0) {obj@data = obj@data[-c(which(is.na(obj@data[i+3]))),]}
    }
  }

  return(obj)})

setMethod('get_model', "Algorithm.Niche.Model", function(obj, ....) {return(obj)})

setMethod('project', "Algorithm.Niche.Model",  function(obj, Env, ...) {
  model = get_model(obj, ...)
  #proj = predict(Env, model) # Previous classic function that don't take correctly factors variable into account
  proj = predict(Env, model,
                 fun = function(model, x){
                   x= as.data.frame(x)
                   for (i in 1:length(Env@layers)) {
                     if(Env[[i]]@data@isfactor) {
                       x[,i] = as.factor(x[,i])
                       levels(x[,i]) = Env[[i]]@data@attributes[[1]]$ID
                     }
                   }
                   return(predict(model, x))
                 })
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = stack(proj)
  return(obj)})

setMethod('evaluate.axes', "Algorithm.Niche.Model", function(obj, thresh = 1001, Env, ...) {
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  for (i in 5:length(obj@data)) {
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    data = obj.axes@data[which(!obj.axes@data$Train),]
    model = get_model(obj.axes,...)
    predicted.values = predict(model, data)
    threshold = optim.thresh(data$Presence, predicted.values, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    evaluation = accuracy(data$Presence, predicted.values, threshold)
    obj@variables.importance[(i-4)] = obj@evaluation$AUC - evaluation$AUC
  }
  if(sum(obj@variables.importance) != 0) {
    obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  } else {
      obj@variables.importance = obj@variables.importance * 100
      }
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

setMethod('sum', 'Algorithm.Niche.Model', function(x, ..., name = NULL, AUCthresh = 0.75, thresh = 1001, format = T, verbose = T, na.rm = F) {
  models = list(x, ...)
  if(format) {for(i in 1:length(models)) {if(!inherits(models[[i]], class(x)[[1]])) {stop('You can only sum models from the same algorithm')}}}
  smodel =  new(class(x)[[1]], 
                projection = stack(reclassify(x@projection[[1]], c(-Inf,Inf,0))),
                data = x@data[1,],
                variables.importance = x@variables.importance)
  smodel@data = smodel@data[-1,]
  smodel@variables.importance[1,] = 0
  
  # Name
  if (!is.null(name)) {name = paste0(name,'.')}
  smodel@name = paste0(name,class(x)[[1]],'.ensemble')
  
  # Datas, Projections, and Variables importance fusion
  sAUC = 0
  kept.model = 0
  for (i in 1:length(models)) {
    if (models[[i]]@evaluation$AUC > AUCthresh) {
      smodel@projection[[1]] = smodel@projection[[1]] + models[[i]]@projection[[1]]*models[[i]]@evaluation$AUC
      smodel@variables.importance = smodel@variables.importance + models[[i]]@variables.importance*models[[i]]@evaluation$AUC
      smodel@data = rbind(smodel@data, models[[i]]@data)
      sAUC = sAUC + models[[i]]@evaluation$AUC
      kept.model = kept.model + 1
    }
  }
  
  # Return NULL if any model is kept
  if( kept.model == 0) {
    if (verbose) {cat('No model were kept with this threshold, Null is return. \n')}
    return(NULL)} else {
      
      smodel@projection[[1]] = smodel@projection[[1]] / sAUC
      names(smodel@projection[[1]]) = 'Probability'
      smodel@variables.importance = smodel@variables.importance / sum(smodel@variables.importance) * 100
      
      # Parameters
      smodel@parameters = x@parameters
      smodel@parameters['kept.model'] = kept.model
      
      # Evaluation
      smodel = evaluate(smodel, thresh = 1001)}})

setMethod('ensemble', 'Algorithm.Niche.Model', function(x, ..., name = NULL, AUCthresh = 0.75, thresh = 1001, uncertainity = T) {
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
    type.model['verbose'] = F
    type.model['name'] = name
    type.model['AUCthresh'] = AUCthresh
    type.model['thresh'] = thresh
   algo.ensemble[type] = do.call(sum, type.model)
  }
  cat('   done. \n')
  
  if (length(algo.ensemble) < 1) {
    cat('No model were kept with this threshold, Null is returned. \n')
    return(NULL)
  } else {
    
    # Sum of algorithm ensemble
    cat('Projcetion, and variables importance computing...')
    algo.list = list()
    for (i in 1:length(algo.ensemble)) {algo.list[[i]] = algo.ensemble[[i]]}
    algo.list['format'] = F
    algo.list['AUCthresh'] = 0
    algo.list['thresh'] = thresh
    sum.algo.ensemble = do.call(sum, algo.list)
    
    # Name
    if (!is.null(name)) {name = paste0(name,'.')}
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
    
    # Binary map
    cat('Binary tranformation...')
    enm@projection[[2]] = reclassify(enm@projection[[1]], c(-Inf,enm@evaluation$threshold,0, enm@evaluation$threshold,Inf,1))
    names(enm@projection[[2]]) = 'Niche'
    cat('   done \n')
    
    # Projections stack
    projections = stack()
    for (i in 1:length(algo.ensemble)) {
      projections = stack(projections, algo.ensemble[[i]]@projection[[1]])
      names(projections[[i]]) = algo.ensemble[[i]]@name
    }
    
    # Algorithms Correlation
    if (uncertainity && length(projections@layers) > 1) {
      cat('Algorithms correlation...')
      enm@algorithm.correlation = as.data.frame(layerStats(projections, 'pearson', na.rm = T)$`pearson correlation coefficient`)
      cat('   done \n')
    }
    
    # Uncertainity map
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
    cat('   done \n')
    
    return(enm)}})

##### GLM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GLM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GLM.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = 1000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GLM.Niche.Model", 
          function(obj, test = 'AIC', epsilon = 1e-08, maxit = 500) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (is.factor(data[,i])) {var = paste0('factor(',var,')')}
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
            }
            model = glm(formula(formula), data = data, test = test, control = glm.control(epsilon = epsilon, maxit = maxit))
            return(model)})

##### GAM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GAM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GAM.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = 1000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GAM.Niche.Model", 
          function(obj, test = 'AIC', epsilon = 1e-08, maxit = 500) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            formula = "Presence ~"
            for (i in 2:length(data)) {
              var = names(data[i])
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
              if (!is.factor(data[,i])) {formula = paste0(formula,' + s(',var,')')}
            }
            model = gam(formula(formula), data = data, test = test, control = gam.control(epsilon = epsilon, maxit = maxit))
            return(model)})

##### MARS Niche Model Class ##### -----

# 1 - Class definition #
setClass('MARS.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "MARS.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = 100
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "MARS.Niche.Model", 
          function(obj, degree = 2) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = earth(Presence ~ ., data = data, degree = 2)
            return(model)})

##### CTA Niche Model Class ##### -----

# 1 - Class definition #
setClass('CTA.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "CTA.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "CTA.Niche.Model", 
          function(obj, final.leave = 1, cv = 3) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = rpart(Presence ~ ., data = data, method = 'class', 
                          control = rpart.control(minbucket = final.leave, xval = cv))
            return(model)})

setMethod('evaluate.axes', "CTA.Niche.Model", function(obj, thresh = 1001, ...) {
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  for (i in 5:length(obj@data)) {
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    data = obj.axes@data[which(!obj.axes@data$Train),]
    model = get_model(obj.axes,...)
    predicted.values = predict(model,data)[,2]
    threshold = optim.thresh(data$Presence, predicted.values, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    evaluation = accuracy(data$Presence, predicted.values, threshold)
    obj@variables.importance[(i-4)] = obj@evaluation$AUC - evaluation$AUC
  }
  obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

##### GBM Niche Model Class ##### -----

# 1 - Class definition #
setClass('GBM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "GBM.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "GBM.Niche.Model", 
          function(obj, trees = 2500, final.leave = 1, cv = 3, thresh.shrink = 1e-03) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = gbm(Presence ~ ., data = data,
                        distribution = 'bernoulli', n.minobsinnode = final.leave,
                        shrinkage = thresh.shrink, bag.fraction = 0.5, 
                        train.fraction = 1, cv.folds = cv, n.trees = trees)
            return(model)})

##### RF Niche Model Class ##### -----

# 1 - Class definition #
setClass('RF.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "RF.Niche.Model", 
          function(obj) {
            PA = list()
#             PA['nb'] = length(obj@data$Presence)
#             PA['strat'] = '2nd'
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "RF.Niche.Model", 
          function(obj, trees = 2500, final.leave = 1) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = randomForest(Presence ~ ., data = data, do.classif = TRUE, ntree = trees, 
                                 nodesize = final.leave, maxnodes = NULL)
            return(model)})

##### MAXENT Niche Model Class ##### -----

# 1 - Class definition #
setClass('MAXENT.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "MAXENT.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = 10000
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "MAXENT.Niche.Model", function(obj, Env) {
  factors = c()
  data = obj@data[which(obj@data$Train),]
  for(i in 4:length(names(obj@data))) {if(is.factor(obj@data[,i])) {factors = c(factors, names(obj@data)[i])}}
  model = maxent(x = Env, p = obj@data[which(data$Presence == 1),1:2], 
                 a = obj@data[which(data$Presence == 0),1:2], factors = factors)
  return(model)})

setMethod('project', "MAXENT.Niche.Model", function(obj, Env, ...) {
  model = get_model(obj, Env)
  proj = predict(Env, model)
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = stack(proj)
  return(obj)})

setMethod('evaluate.axes', "MAXENT.Niche.Model", function(obj, thresh = 1001, Env, ...) {
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  for (i in 5:length(obj@data)) {
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    model = get_model(obj.axes, Env[[-(i-4)]])
    obj.axes@data = obj.axes@data[which(!obj.axes@data$Train),]
    predicted.values = predict(model, obj.axes@data)
    threshold = optim.thresh(obj.axes@data$Presence, predicted.values, thresh)
    threshold = mean(threshold$`max.sensitivity+specificity`)
    evaluation = accuracy(obj.axes@data$Presence, predicted.values, threshold)
    obj@variables.importance[(i-4)] = obj@evaluation$AUC - evaluation$AUC
  }
  obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

##### ANN Niche Model Class ##### -----

# 1 - Class definition #
setClass('ANN.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "ANN.Niche.Model", 
          function(obj) {
            PA = list()
            PA['nb'] = length(obj@data$Presence)
            PA['strat'] = 'random'
            return(PA)})

setMethod('get_model', "ANN.Niche.Model", function(obj, maxit = 500) {
            data = obj@data[which(obj@data$Train),]
            data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
            model = nnet(Presence ~ ., data = data, size = 6, maxit = maxit)
            return(model)})

##### SVM Niche Model Class ##### -----

# 1 - Class definition #
setClass('SVM.Niche.Model',
         contains = 'Algorithm.Niche.Model')

# 2 - Methods definition #
setMethod('get_PA', "SVM.Niche.Model", function(obj) {
  PA = list()
  PA['nb'] = length(obj@data$Presence)
  PA['strat'] = 'random'
  return(PA)})

setMethod('get_model', "SVM.Niche.Model", function(obj, epsilon = 1e-08, cv = 3) {
  data = obj@data[which(obj@data$Train),]
  data = data[-c(which(names(data) == 'X'),which(names(data) == 'Y'),which(names(data) == 'Train'))]
  model = svm(Presence ~ ., data = data, type = 'eps-regression', 
              gamma = 1/(length(data)-1), kernel = 'radial', epsilon = epsilon, cross = cv)
  return(model)})

##### Ensemble Niche Model Class ##### -----

# 1 - Class definition #
setClass('Ensemble.Niche.Model',
         contains = 'Niche.Model',
         representation(uncertainity = 'Raster',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame'),
         prototype(uncertainity = raster(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame()))

# 2 - Class creation function #
Ensemble.Niche.Model <- function(name = character(), 
                                 projection = stack(), 
                                 evaluation = data.frame(), 
                                 variables.importance = data.frame(),
                                 data = data.frame(),
                                 parameters = list(),
                                 uncertainity = raster(),
                                 algorithm.correlation = data.frame(),
                                 algorithm.evaluation = data.frame()) {
  return(new('Ensemble.Niche.Model', 
             name = name, 
             projection = projection, 
             evaluation = evaluation, 
             variables.importance = variables.importance, 
             data = data, 
             parameters = parameters,
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation))}

load.enm = function (directory = getwd()) {
  a = try(read.csv(paste0(directory,'/Tables/AlgoCorr'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty !')
    a = data.frame()
    }
  enm = Ensemble.Niche.Model(name = as.character(read.csv(paste0(directory,'/Tables/Name'))[1,2]),
                             projection = stack(raster(paste0(directory,'/Rasters/Probability.tif')), raster(paste0(directory,'/Rasters/Niche.tif'))),
                             uncertainity = try(raster(paste0(directory,'/Rasters/Uncertainity.tif'))),
                             evaluation = read.csv(paste0(directory,'/Tables/ENMeval'), row.names = 1),
                             algorithm.evaluation  = read.csv(paste0(directory,'/Tables/AlgoEval'), row.names = 1),
                             algorithm.correlation = a,
                             data = read.csv(paste0(directory,'/Tables/Data'), row.names = 1),
                             variables.importance = read.csv(paste0(directory,'/Tables/VarImp'), row.names = 1))
  return(enm)
}

# 3 - Methods definition #
setMethod('save.enm', 'Ensemble.Niche.Model', function (enm, directory = getwd()) {
            
            cat('Saving ensemble model results \n')
            # Directories creation
            dir.create(path = paste0(directory, "/Rasters"))
            dir.create(path = paste0(directory, "/Tables"))
            
            # Raster saving
            cat('   rasters ...')
            writeRaster(enm@projection[[1]], paste0(directory,'/Rasters/Probability'), 'GTiff', overwrite = T)
            writeRaster(enm@projection[[2]], paste0(directory,'/Rasters/Niche'), 'GTiff', overwrite = T)
            writeRaster(enm@uncertainity, paste0(directory,'/Rasters/Uncertainity'), 'GTiff', overwrite = T)
            cat('saved \n')
            
            # Tables saving
            cat('   tables ...')
            write.csv(enm@evaluation, paste0(directory,'/Tables/ENMeval'))
            write.csv(enm@algorithm.evaluation, paste0(directory,'/Tables/AlgoEval'))
            write.csv(enm@algorithm.correlation, paste0(directory,'/Tables/AlgoCorr'))
            write.csv(enm@variables.importance, paste0(directory,'/Tables/VarImp'))
            write.csv(enm@data, paste0(directory,'/Tables/Data'))
            write.csv(enm@name, paste0(directory,'/Tables/Name'))
            cat('saved \n \n')
          })

setMethod('stacking', 'Ensemble.Niche.Model', function(enm, ..., name = NULL) {
  enms = list(enm, ...)
  if (length(enms) < 2) {stop('You neeed more than one enm to do stackings')}
  cat('Stack creation... \n')
  stack = Stack.Species.Ensemble.Niche.Model(diversity.map = reclassify(enm@projection[[1]], c(-Inf,Inf,0)), 
                                             uncertainity = reclassify(enm@uncertainity, c(-Inf,Inf,NA)))
  
  # Name
  cat('   naming...')
  if (!is.null(name)) {name = paste0(name,'.')}
  stack@name = paste0(name,'Stack.Species.Ensemble.Niche.Model')
  cat(' done. \n')
  
  # Diversity map
  cat('   diversity mapping...')
  for (i in 1:length(enms)) {stack@diversity.map = stack@diversity.map + enms[[i]]@projection[[1]]}
  names(stack@diversity.map) = 'diversity'
  cat(' done. \n')
  
  # Uncertainity map
  cat('   uncertainity mapping...')
  uncertainities = stack()
  for (i in 1:length(enms)) {uncertainities = stack(uncertainities, enms[[i]]@uncertainity)}
  stack@uncertainity = calc(uncertainities, mean)
  names(stack@uncertainity) = 'uncertainity'
  cat(' done.')
  
  # Evaluation
  cat('   evaluating...')
  stack@evaluation = enm@evaluation
  for (i in 2:length(enms)) {stack@evaluation = rbind(stack@evaluation, enms[[i]]@evaluation)}
  a = stack@evaluation[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in 1:length(stack@evaluation)) {a[i] = c(mean(stack@evaluation[,i]), var(stack@evaluation[i]))}
  stack@evaluation = a
  cat(' done. \n')
  
  # Variables Importance
  cat('   comparing variables importnace...')
  stack@variables.importance = enm@variables.importance
  for (i in 2:length(enms)) {stack@variables.importance = rbind(stack@variables.importance, enms[[i]]@variables.importance)}
  a = stack@variables.importance[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in 1:length(stack@variables.importance)) {a[i] = c(mean(stack@variables.importance[,i]), var(stack@variables.importance[i]))}
  stack@variables.importance = a
  cat(' done. \n')
  
  # Algorithm Correlation
  cat('   comparing algorithms correlation...')
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
  
  stack@enms = enms
  
  return(stack)})

##### Stack Species Ensemble Niche Model Class ##### -----

# 1 - Class definition #
setClass('Stack.Species.Ensemble.Niche.Model',
         representation(name = 'character',
                        diversity.map = 'Raster',
                        uncertainity = 'Raster',
                        evaluation = 'data.frame',
                        variables.importance = 'data.frame',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame',
                        enms = 'list'),
         prototype(name = character(),
                   diversity.map = raster(),
                   uncertainity = raster(),
                   evaluation = data.frame(),
                   variables.importance = data.frame(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame(),
                   enms = list()))

# 2 - Class creation function #
Stack.Species.Ensemble.Niche.Model <- function(name = character(),
                                               diversity.map = raster(),
                                               uncertainity = raster(),
                                               evaluation = data.frame(),
                                               variables.importance = data.frame(),
                                               algorithm.correlation = data.frame(),
                                               algorithm.evaluation = data.frame(),
                                               enms = list()) {
  return(new('Stack.Species.Ensemble.Niche.Model', 
             name = name, 
             diversity.map = diversity.map, 
             evaluation = evaluation, 
             variables.importance = variables.importance, 
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             enms = enms))}

load.sp.enms = function (directory = getwd()) {
  enm = Ensemble.Niche.Model(name = as.character(read.csv(paste0(directory,'/Tables/Name'))[1,2]),
                             projection = stack(raster(paste0(directory,'/Rasters/Probability.tif')), raster(paste0(directory,'/Rasters/Niche.tif'))),
                             uncertainity = raster(paste0(directory,'/Rasters/Uncertainity.tif')),
                             evaluation = read.csv(paste0(directory,'/Tables/ENMeval')),
                             algorithm.evaluation  = read.csv(paste0(directory,'/Tables/AlgoEval')),
                             algorithm.correlation = read.csv(paste0(directory,'/Tables/AlgoCorr')),
                             data = read.csv(paste0(directory,'/Tables/Data')),
                             variables.importance = read.csv(paste0(directory,'/Tables/VarImp')))
  return(enm)
}

# 3 - Methods definition #
setMethod('save.enm', 'Stack.Species.Ensemble.Niche.Model', function (enm, directory = getwd()) {
  
  cat('Saving ensemble model results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/Rasters"))
  dir.create(path = paste0(directory, "/Tables"))
  
  # Raster saving
  cat('   rasters ...')
  writeRaster(enm@projection[[1]], paste0(directory,'/Rasters/Probability'), 'GTiff', overwrite = T)
  writeRaster(enm@projection[[2]], paste0(directory,'/Rasters/Niche'), 'GTiff', overwrite = T)
  writeRaster(enm@uncertainity, paste0(directory,'/Rasters/Uncertainity'), 'GTiff', overwrite = T)
  cat('saved \n')
  
  # Tables saving
  cat('   tables ...')
  write.csv(enm@evaluation, paste0(directory,'/Tables/ENMeval'))
  write.csv(enm@algorithm.evaluation, paste0(directory,'/Tables/AlgoEval'))
  write.csv(enm@algorithm.correlation, paste0(directory,'/Tables/AlgoCorr'))
  write.csv(enm@variables.importance, paste0(directory,'/Tables/VarImp'))
  write.csv(enm@data, paste0(directory,'/Tables/Data'))
  write.csv(enm@name, paste0(directory,'/Tables/Name'))
  cat('saved \n \n')
})

setMethod('plot', 'Stack.Species.Ensemble.Niche.Model', function(x, y, ...) {
  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Stacked species", tabName = "stack", icon = icon("dashboard")),
        menuItemOutput('ENMS')
      )
    ),
    dashboardBody(
      tabItems(
        # Main page beginning #
        tabItem(tabName = "stack",
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                                   plotOutput('Diversity', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)), title = 'Diversity'),
                         tabPanel(plotOutput('Uncertainity'), title = 'Uncertainity')
                  ),
                  tabBox(title = 'Variables importance',
                         tabPanel(plotOutput('varimp.barplot'), title = 'Barplot'), 
                         tabPanel(tableOutput('varimp.table'), title = 'Table')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('evaluation.barplot'), title = 'Barplot'), 
                         tabPanel(tableOutput('evaluation.table'), title = 'Table')
                  ),
                  if (length(x@algorithm.correlation) > 0) {
                    tabBox(title = 'Algorithms correlation',
                           tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'), 
                           tabPanel(tableOutput('algo.corr.table'), title = 'Table')
                    )
                  }
                )
        )#,
#         # Main page end #
#         # ENM page beginning #
#         for (i in 1:length(x@enms)-1) {
#           tabItem(tabName = ),
#                   h2('enm'))
#         },
#         tabItem(tabName = 'enm2',
#                 h2('enm2')
#         )
      )
    )
  )
  
  server <- function(input, output) {
    
    # Main page beginning #
    if (length(x@algorithm.correlation) > 0) {
      x@algorithm.correlation[upper.tri(x@algorithm.correlation)] = NA
      names(x@algorithm.correlation) = row.names(x@algorithm.correlation)
      x@algorithm.correlation = x@algorithm.correlation[-length(x@algorithm.correlation)]
      x@algorithm.correlation = x@algorithm.correlation[-length(x@algorithm.correlation)]
    }
    # Maps
    # Single zoomable plot
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$unzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    output$Diversity <- renderPlot({
      if (!is.null(ranges$x)) {diversity = crop(x@diversity.map, c(ranges$x, ranges$y))} else {diversity = x@diversity.map}
      plot(diversity, 
           main = paste('AUC :',round(x@evaluation$AUC[1],3),'  Kappa',round(x@evaluation$Kappa[1],3)),
           xlab = 'Longitude (°)',
           ylab = 'Latitude (°)')
    })
    output$uncertainity <- renderPlot({
      if (!is.null(ranges$x)) {uncert = crop(x@duncertainity, c(ranges$x, ranges$y))} else {uncert = x@uncertainity}
      plot(uncert, main = paste('AUC :',round(x@evaluation$AUC[1],3),'  Kappa',round(x@evaluation$Kappa[1],3)))})
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation = x@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$evaluation.table <- renderTable({x@algorithm.evaluation[c(2,4:8)]})
    if (length(x@algorithm.correlation) > 0) {
      # Algorithms correlation
      output$algo.corr.table <- renderTable({x@algorithm.correlation})
      output$algo.corr.heatmap <- renderPlot({
        m <- as.matrix(x@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,3), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
      })
    }
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@variables.importance[-1]))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = row.names(varimp), las = 2)
    })
    output$varimp.table <- renderTable({x@variables.importance[-1]})
    
    output$ENMS <- renderMenu({
      enms = list()
      for (i in 1:length(x@enms)) {enms[[i]] = menuSubItem(x@enms[[i]]@name, tabName = paste0('enm',i))}
      sidebarMenu(.list = enms)
    })
  }
    
  shinyApp(ui, server)
  })
