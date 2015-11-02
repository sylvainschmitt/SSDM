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
setGeneric('evaluate', function(obj, thresh = 1001, metric = 'AUC') {return(standardGeneric('evaluate'))})
setGeneric('get_PA', function(obj) {return(standardGeneric('get_PA'))})
setGeneric('PA.select', function(obj, Env, ...) {return(standardGeneric('PA.select'))})
setGeneric('data.values', function(obj, Env, na.rm = T) {return(standardGeneric('data.values'))})
setGeneric('get_model', function(obj, ...) {return(standardGeneric('get_model'))})
setGeneric('project', function(obj, Env, ...) {return(standardGeneric('project'))})
setGeneric('evaluate.axes', function(obj, thresh = 1001, Env, axes.metric = 'AUC', ...) {return(standardGeneric('evaluate.axes'))})
setGeneric('ensemble', function(x, ..., name = NULL,ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = T, thresh = 1001, uncertainity = T) {return(standardGeneric('ensemble'))})
setGeneric('save.enm', function (enm, ...) {return(standardGeneric('save.enm'))})
setGeneric('save.stack', function (stack, ...) {return(standardGeneric('save.stack'))})
setGeneric('stacking', function(enm, ...) {return(standardGeneric('stacking'))})

##### Niche Model Class ##### -----

# 1 - Class definition #
setClass('Niche.Model',
         representation(name = 'character', 
                        projection = 'Raster', 
                        evaluation = 'data.frame', 
                        variables.importance = 'data.frame',
                        data = 'data.frame',
                        parameters = 'data.frame'),
         prototype(name = character(), 
                   projection = raster(), 
                   evaluation = data.frame(), 
                   variables.importance = data.frame(),
                   data = data.frame(),
                   parameters = data.frame()))

# 2 - Methods definition #
setMethod("evaluate", "Niche.Model", function(obj, thresh = 1001, metric = 'SES') {
  obj@parameters$metric = metric
  data = obj@data[which(!obj@data$Train),]
  predicted.values = extract(obj@projection, data[c('X','Y')])
  predicted.values[which(is.na(predicted.values))] = 0
  # Threshold computing
  metric = switch(metric,
                  'Kappa' = 'maxKappa',
                  'CCR' = 'max.prop.correct',
                  'TSS' = 'max.sensitivity+specificity',
                  'SES' = 'sensitivity=specificity',
                  'LW' = 'min.occurence.prediction',
                  'ROC' = 'min.ROC.plot.distance')
  threshold = optim.thresh(data$Presence, predicted.values, thresh)
  threshold = mean(threshold[[which(names(threshold) == metric)]])
  
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
        if (!is.null(ranges$x)) {proba.map = crop(x@projection, c(ranges$x, ranges$y))} else {proba.map = x@projection}
        plot(proba.map, 
             main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)),
             xlab = 'Longitude (°)',
             ylab = 'Latitude (°)')
        points(x@data$X[which(x@data$Presence == 1)], 
               x@data$Y[which(x@data$Presence == 1)], 
               pch = 16, cex = 0.7)
      })
      output$niche <- renderPlot({
        if (!is.null(ranges$x)) {
          niche.map = crop(reclassify(x@projection, c(-Inf,x@evaluation$threshold,0, x@evaluation$threshold,Inf,1)), c(ranges$x, ranges$y))
        } else {niche.map = reclassify(x@projection, c(-Inf,x@evaluation$threshold,0, x@evaluation$threshold,Inf,1))}
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
                                  projection = raster(), 
                                  evaluation = data.frame(), 
                                  variables.importance = data.frame(),
                                  data = data.frame(),
                                  parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  object.class = paste0(algorithm,'.Niche.Model')
  return(new(object.class, name = name, projection = projection, evaluation = evaluation, variables.importance = variables.importance, data = data, parameters = parameters))
}

# 3 - Methods definition #
setMethod('get_PA', "Algorithm.Niche.Model", function(obj) {return(obj)})

setMethod('PA.select', "Algorithm.Niche.Model", function(obj, Env, PA = NULL, train.frac = 0.7) {
  if (is.null(PA)) {
    PA = get_PA(obj)
    obj@parameters$PA = 'default'
  } else {
    obj@parameters$PA = paste0(as.character(PA$nb),'.',as.character(PA$strat))
  }
  
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
                       x[,i] = droplevels(x[,i])
                       levels(x[,i]) = Env[[i]]@data@attributes[[1]]$ID
                     }
                   }
                   return(predict(model, x))
                 })
  # Rescaling projection
  proj = reclassify(proj, c(-Inf,0,0))
  proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  return(obj)})

setMethod('evaluate.axes', "Algorithm.Niche.Model", function(obj, thresh = 1001, Env, 
                                                             axes.metric = 'AUC', ...) {
  obj@parameters$axes.metric = axes.metric
  obj@variables.importance = data.frame(matrix(nrow = 1, ncol = (length(obj@data)-4)))
  names(obj@variables.importance) = names(obj@data)[5:length(obj@data)]
  if (axes.metric == 'Pearson') {
    data = obj@data[which(!obj@data$Train),]
    o.predicted.values = extract(obj@projection, data[c('X','Y')]) # original model predicted values
  }

  for (i in 5:length(obj@data)) {
    # Get model predictions without one axis reeated for all axis
    obj.axes = obj
    obj.axes@data = obj.axes@data[-i]
    data = obj.axes@data[which(!obj.axes@data$Train),]
    model = get_model(obj.axes,...)
    predicted.values = predict(model, data)
    if (axes.metric != 'Pearson') {
      threshold = optim.thresh(data$Presence, predicted.values, thresh)
      threshold = mean(threshold$`max.sensitivity+specificity`)
      evaluation = accuracy(data$Presence, predicted.values, threshold)
      obj@variables.importance[1,(i-4)] = obj@evaluation[1,which(names(obj@evaluation) == axes.metric)]  - evaluation[1,which(names(evaluation) == axes.metric)]
    } else {
      obj@variables.importance[(i-4)] = cor(predicted.values, o.predicted.values)
    }
  }
  
  # Variable importance normalization (%)
  if(sum(obj@variables.importance) == 0) {
    all.null = T
    for (i in 1:length(obj@variables.importance[1,])) {if(obj@variables.importance[1,i] != 0) {all.null = F}}
    if (all.null) {
      obj@variables.importance[1,] = 100 / length(obj@variables.importance)
    } else {
      obj@variables.importance = obj@variables.importance * 100
    }
  } else {
    obj@variables.importance = obj@variables.importance / sum(obj@variables.importance) * 100
  }
  row.names(obj@variables.importance) = "Axes.evaluation"
  return(obj)})

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
      smodel = evaluate(smodel, thresh = 1001)}})

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
      smodel@parameters$ensemble.metric = text.ensemble.metric 
      smodel@parameters$ensemble.thresh = text.ensemble.thresh
      smodel@parameters$weight = weight
      smodel@parameters$kept.model = kept.model
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
              if (i != 2) {formula = paste(formula,'+',var)} else {formula = paste(formula,var)}
            }
            model = glm(formula(formula), data = data, test = test, control = glm.control(epsilon = epsilon, maxit = maxit))
            for(i in 1:length(data)) {
              if(is.factor(data[,i])) {
                model$xlevels[[which(names(model$xlevels) == paste0(names(data)[i]))]] = levels(data[,i])
              }
            }
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
#             for(i in 1:length(data)) {
#               if(is.factor(data[,i])) {
#                 model$xlevels[[which(names(model$xlevels) == names(data)[i])]] = levels(data[,i])
#               }
#             }
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
  obj@projection = proj
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
                                 projection = raster(), 
                                 evaluation = data.frame(), 
                                 variables.importance = data.frame(),
                                 data = data.frame(),
                                 uncertainity = raster(),
                                 algorithm.correlation = data.frame(),
                                 algorithm.evaluation = data.frame(),
                                 parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Ensemble.Niche.Model', 
             name = name, 
             projection = projection, 
             evaluation = evaluation, 
             variables.importance = variables.importance, 
             data = data, 
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             parameters = parameters))}

load.enm = function (name, directory = getwd()) {
  directory = paste0(directory, '/', name)
  a = try(read.csv(paste0(directory,'/Tables/AlgoCorr'), row.names = 1))
  if (inherits(a, 'try-error')) {
    cat('Algorithm correlation table empty !')
    a = data.frame()
    }
  enm = Ensemble.Niche.Model(name = as.character(read.csv(paste0(directory,'/Tables/Name'))[1,2]),
                             projection = raster(paste0(directory,'/Rasters/Probability.tif')),
                             uncertainity = try(raster(paste0(directory,'/Rasters/Uncertainity.tif'))),
                             evaluation = read.csv(paste0(directory,'/Tables/ENMeval'), row.names = 1),
                             algorithm.evaluation  = read.csv(paste0(directory,'/Tables/AlgoEval'), row.names = 1),
                             algorithm.correlation = a,
                             data = read.csv(paste0(directory,'/Tables/Data'), row.names = 1),
                             variables.importance = read.csv(paste0(directory,'/Tables/VarImp'), row.names = 1),
                             parameters = read.csv(paste0(directory,'/Tables/Parameters'), row.names = 1))
  return(enm)
}

# 3 - Methods definition #
setMethod('save.enm', 'Ensemble.Niche.Model', function (enm, 
                                                        name = strsplit(enm@name, '.', fixed = T)[[1]][1], 
                                                        directory = getwd()) {
  
  cat('Saving ensemble model results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/", name))
  dir.create(path = paste0(directory, "/", name,"/Rasters"))
  dir.create(path = paste0(directory, "/", name, "/Tables"))
  
  # Raster saving
  cat('   rasters ...')
  writeRaster(enm@projection[[1]], paste0(directory, "/", name, '/Rasters/Probability'), 'GTiff', overwrite = T)
  writeRaster(enm@uncertainity, paste0(directory, "/", name, '/Rasters/Uncertainity'), 'GTiff', overwrite = T)
  cat('saved \n')
  
  # Tables saving
  cat('   tables ...')
  write.csv(enm@evaluation, paste0(directory, "/", name, '/Tables/ENMeval'))
  write.csv(enm@algorithm.evaluation, paste0(directory, "/", name, '/Tables/AlgoEval'))
  write.csv(enm@algorithm.correlation, paste0(directory, "/", name, '/Tables/AlgoCorr'))
  write.csv(enm@variables.importance, paste0(directory, "/", name, '/Tables/VarImp'))
  write.csv(enm@data, paste0(directory, "/", name, '/Tables/Data'))
  write.csv(enm@name, paste0(directory, "/", name, '/Tables/Name'))
  write.csv(enm@parameters, paste0(directory, "/", name, '/Tables/Parameters'))
  cat('saved \n \n')
})

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
                        enms = 'list',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   diversity.map = raster(),
                   uncertainity = raster(),
                   evaluation = data.frame(),
                   variables.importance = data.frame(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame(),
                   enms = list(),
                   parameters = data.frame()))

# 2 - Class creation function #
Stack.Species.Ensemble.Niche.Model <- function(name = character(),
                                               diversity.map = raster(),
                                               uncertainity = raster(),
                                               evaluation = data.frame(),
                                               variables.importance = data.frame(),
                                               algorithm.correlation = data.frame(),
                                               algorithm.evaluation = data.frame(),
                                               enms = list(),
                                               parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Stack.Species.Ensemble.Niche.Model', 
             name = name, 
             diversity.map = diversity.map, 
             evaluation = evaluation, 
             variables.importance = variables.importance, 
             uncertainity = uncertainity,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             enms = enms,
             parameters = parameters))}

load.stack = function (name = 'Stack', directory = getwd()) {
  directory = paste0(directory, '/', name)
  stack = Stack.Species.Ensemble.Niche.Model(name = as.character(read.csv(paste0(directory,'/Results/Tables/Name'))[1,2]),
                                             diversity.map = raster(paste0(directory,'/Results/Rasters/Diversity.tif')),
                                             uncertainity = raster(paste0(directory,'/Results/Rasters/Uncertainity.tif')),
                                             evaluation = read.csv(paste0(directory,'/Results/Tables/StackEval'), row.names = 1),
                                             variables.importance = read.csv(paste0(directory,'/Results/Tables/VarImp'), row.names = 1),
                                             algorithm.correlation = read.csv(paste0(directory,'/Results/Tables/AlgoCorr'), row.names = 1),
                                             algorithm.evaluation = read.csv(paste0(directory,'/Results/Tables/AlgoEval'), row.names = 1),
                                             enms = list(),
                                             parameters = read.csv(paste0(directory,'/Results/Tables/Parameters'), row.names = 1))
  enms = list.dirs(directory, recursive = F, full.names = F)
  enms = enms[-which(enms == 'Results')]
  for (i in 1:length(enms)) {stack@enms[[i]] = load.enm(enms[i], directory = directory)}
  return(stack)
}

# 3 - Methods definition #
setMethod('save.stack', 'Stack.Species.Ensemble.Niche.Model', function (stack, name = 'Stack', directory = getwd()) {
  
  cat('Saving stack species model results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/", name))
  directory = paste0(directory, "/", name)
  dir.create(path = paste0(directory, "/", 'Results'))
  dir.create(path = paste0(directory, "/", 'Results',"/Rasters"))
  dir.create(path = paste0(directory, "/", 'Results', "/Tables"))
  
  # Raster saving
  cat('   rasters ...')
  writeRaster(stack@diversity.map, paste0(directory, "/", 'Results', '/Rasters/Diversity'), 'GTiff', overwrite = T)
  writeRaster(stack@uncertainity, paste0(directory, "/", 'Results', '/Rasters/Uncertainity'), 'GTiff', overwrite = T)
  cat('saved \n')
  
  # Tables saving
  cat('   tables ...')
  write.csv(stack@evaluation, paste0(directory, "/", 'Results', '/Tables/StackEval'))
  write.csv(stack@algorithm.evaluation, paste0(directory, "/", 'Results', '/Tables/AlgoEval'))
  write.csv(stack@algorithm.correlation, paste0(directory, "/", 'Results', '/Tables/AlgoCorr'))
  write.csv(stack@variables.importance, paste0(directory, "/", 'Results', '/Tables/VarImp'))
  write.csv(stack@name, paste0(directory, "/", 'Results', '/Tables/Name'))
  write.csv(stack@parameters, paste0(directory, "/", 'Results', '/Tables/Parameters'))
  cat('saved \n\n')
  
  # ENMS saving
  cat('   enms ... \n\n')
  for (i in 1:length(stack@enms)) {save.enm(stack@enms[[i]], directory = directory)}
  cat('saved \n \n')
})

setMethod('plot', 'Stack.Species.Ensemble.Niche.Model', function(x, y, ...) {
  choices = list()
  for (i in 1:length(x@enms)) {choices[[i]] = strsplit(x@enms[[i]]@name, '.', fixed = T)[[1]][1]}
  if (inherits(x@enms[[1]], 'Algorithm.Niche.Model')) {full = F} else {full = T}
  
  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Stacked species", tabName = "stack", icon = icon("dashboard")),
        menuItem("Ensemble model", tabName = "enm", icon = icon("pagelines")),
        selectInput('enmchoice', 'Ensemble model specie :', choices, selected = NULL, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      )
    ),
    dashboardBody(
      tabItems(
        # Main page beginning #
        tabItem(tabName = "stack",
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                                   plotOutput('Diversity', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)),
                                   textOutput('diversity.info'),
                                   title = 'Diversity'),
                         tabPanel(plotOutput('Uncertainity'), title = 'Uncertainty'),
                         tabPanel(tableOutput('summary'), title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('varimp.barplot'), 
                                  textOutput('varimp.info'),
                                  title = 'Barplot'), 
                         tabPanel(tableOutput('varimp.table'), title = 'Table')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('evaluation.barplot'), 
                                  textOutput('evaluation.info'),
                                  title = 'Barplot'), 
                         tabPanel(tableOutput('evaluation.table'), title = 'Table')
                  ),
                  if (length(x@algorithm.correlation) > 0) {
                    tabBox(title = 'Algorithm correlation',
                           tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'), 
                           tabPanel(tableOutput('algo.corr.table'), title = 'Table')
                    )
                  }
                )
        ),
        # Main page end #
        # ENM page beginning #
        tabItem(tabName = 'enm',
          fluidRow(
            tabBox(title = 'Maps',
                   tabPanel( actionButton('enmunzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                             plotOutput('probability', dblclick = "proba_dblclick", brush = brushOpts(id = "proba_brush", resetOnNew = TRUE)),
                             title = 'Habitat suitability'),
                   tabPanel(plotOutput('niche'), 
                            textOutput('enm.binary.info'),
                            title = 'Binary map'),
                   tabPanel(plotOutput('enm.uncertainity'), title = 'Uncertainty'),
                   tabPanel(tableOutput('enm.summary'), title = 'Summary')
            ),
            tabBox(title = 'Variable importance',
                   tabPanel(plotOutput('enm.varimp.barplot'), 
                            textOutput('enm.varimp.info'),
                            title = 'Barplot'), 
                   tabPanel(tableOutput('enm.varimp.table'), title = 'Table')
            )
          ),
          fluidRow(
            tabBox(title = 'Model evaluation',
                   tabPanel(plotOutput('enm.evaluation.barplot'), 
                            textOutput('enm.evaluation.info'),
                            title = 'Barplot'), 
                   tabPanel(tableOutput('enm.evaluation.table'), title = 'Table')
            ),
            
            tabBox(title = 'Algorithm correlation',
                   tabPanel(plotOutput('enm.algo.corr.heatmap'), title = 'Heatmap'), 
                   tabPanel(tableOutput('enm.algo.corr.table'), title = 'Table')
            )
          )
        )
        # ENM page ending #
      )
    )
  )
  
  server <- function(input, output) {
    
    # Main page beginning #
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
    eval = character()
    ensemble.metric = strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
    for (i in 1:length(ensemble.metric)) {
      eval = paste(eval, paste(ensemble.metric[i],':',round(x@evaluation[1,which(names(x@evaluation) == ensemble.metric[i])], digits = 3)))
      if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
    }
    output$Diversity <- renderPlot({
      if (!is.null(ranges$x)) {diversity = crop(x@diversity.map, c(ranges$x, ranges$y))} else {diversity = x@diversity.map}
      plot(diversity, 
           main = eval,
           xlab = 'Longitude (°)',
           ylab = 'Latitude (°)',
           legend.args=list(text='Local \nspecies \nrichness', font = 3, line = 1))
    })
    output$Uncertainity <- renderPlot({
      if (!is.null(ranges$x)) {uncert = crop(x@uncertainity, c(ranges$x, ranges$y))} else {uncert = x@uncertainity}
      plot(uncert, main = eval, legend.args=list(text='Models \nvariance', font = 3, line = 1))})
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation = x@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / x@parameters$rep
      metrics = '% kept.model'
      metrics.nb = c(which(names(evaluation) == 'kept.model'))
      for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        metrics = c(metrics, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i])
        metrics.nb = c(metrics.nb, which(names(evaluation) == strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i]))
      }
      table <- t(evaluation[metrics.nb])
      barplot(table, col = rainbow(length(metrics)), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', metrics, fill = rainbow(length(metrics)))
    })
    output$evaluation.table <- renderTable({x@algorithm.evaluation[c(2,4:8)]})
    if (length(x@algorithm.correlation) > 0) {
      # Algorithms correlation
      output$algo.corr.table <- renderTable({x@algorithm.correlation})
      output$algo.corr.heatmap <- renderPlot({
        m <- as.matrix(x@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,2), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
      })
    }
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@variables.importance))
      names(varimp) = 'Axes.evaluation'
      bar = barplot(varimp$Axes.evaluation, names.arg = strwrap(row.names(varimp)), 
                    ylim = c(0,(max(varimp$Axes.evaluation)+max((varimp[2])))),
                    las = 2, ylab = 'Variable relative contribution (%)')
      arrows(bar,varimp$Axes.evaluation+varimp[,2], bar, varimp$Axes.evaluation-varimp[,2], angle=90, code=3, length=0.1)
    })
    output$varimp.table <- renderTable({x@variables.importance})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 5, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurences type', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
      # 'Original number of species'
      algo.info = character()
      for (i in 1:length(strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (x@parameters$PA) {PA = 'default'}
      summary$Summary = c(x@parameters$data, length(x@enms), algo.info, x@parameters$rep, PA)
      if(!is.null(x@parameters$sp.nb.origin)) {
        summary = rbind(summary, 
                        data.frame(Summary = x@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$diversity.info <- renderText({
      x@parameters$method = switch(x@parameters$method,
                                   'P' = 'summing habitat suitability map probabilities.',
                                   'T' = paste('summing habitat suitability binary map after thresholding with',x@parameters$metric),
                                   'B' = paste('drawing repeatdly bernoulli repetitions with',x@parameters$rep.B))
      text = paste('Local species richness map realized by', x@parameters$method)
      text
    })
    output$varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(x@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, x@parameters$axes.metric[i])
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', x@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', x@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (x@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # Main page ending #
    
    # ENM beginning #
    # Maps
    # Single zoomable plot
    observeEvent(input$proba_dblclick, {
      brush <- input$proba_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$enmunzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    
    output$probability <- renderPlot({
      if (!is.null(ranges$x)) {proba = crop(x@enms[[which(choices == input$enmchoice)]]@projection, c(ranges$x, ranges$y))} else {proba = x@enms[[which(choices == input$enmchoice)]]@projection}
      plot(proba, 
           main = paste('AUC :',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$AUC,3),'  Kappa',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$Kappa,3)),
           xlab = 'Longitude (°)',
           ylab = 'Latitude (°)',
           legend.args=list(text='Presence\nprobability', font = 3, line = 1))
      points(x@enms[[which(choices == input$enmchoice)]]@data$X[which(x@enms[[which(choices == input$enmchoice)]]@data$Presence == 1)], 
             x@enms[[which(choices == input$enmchoice)]]@data$Y[which(x@enms[[which(choices == input$enmchoice)]]@data$Presence == 1)], 
             pch = 16, cex = 0.7)
    })
    output$niche <- renderPlot({
      niche.map = reclassify(x@enms[[which(choices == input$enmchoice)]]@projection, c(-Inf,x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold,0, x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold,Inf,1))
      if (!is.null(ranges$x)) {niche.map = crop(niche.map, c(ranges$x, ranges$y))}
      plot(niche.map, main = paste('AUC :',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$AUC,3),'  Kappa',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$Kappa,3)))})
    output$enm.uncertainity <- renderPlot({
      if (!is.null(ranges$x)) {uncert.map = crop(x@enms[[which(choices == input$enmchoice)]]@uncertainity, c(ranges$x, ranges$y))} else {uncert.map = x@enms[[which(choices == input$enmchoice)]]@uncertainity}
      plot(uncert.map, main = paste('AUC :',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$AUC,3),'  Kappa',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$Kappa,3)), legend.args=list(text='Models \nvariance', font = 3, line = 1))})
    # Evaluation
    output$enm.evaluation.barplot <- renderPlot({
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.', fixe = T)[[1]][2]}
      evaluation = x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$enm.evaluation.table <- renderTable({
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.Niche')[[1]][1]}
      x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation[c(2,4:8)]
    })
    # Algorithms correlation
    output$enm.algo.corr.table <- renderTable({
      if (length(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation) > 0) {
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation[upper.tri(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation, diag = T)] = NA
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation
      }
    })
    output$enm.algo.corr.heatmap <- renderPlot({
      if (length(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation) > 0) {
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation[upper.tri(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation, diag = T)] = NA
        m <- as.matrix(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,3), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
      }
    })
    # Variable importance
    output$enm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@enms[[which(choices == input$enmchoice)]]@variables.importance[-1]))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = row.names(varimp), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$enm.varimp.table <- renderTable({x@enms[[which(choices == input$enmchoice)]]@variables.importance[-1]})
    # Parameters
    output$enm.summary <- renderTable({
      summary = data.frame(matrix(nrow = 5, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurences type', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
      algo.info = character()
      for (i in 1:length(strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (x@enms[[which(choices == input$enmchoice)]]@parameters$PA) {PA = 'default'}
      summary$Summary = c(x@enms[[which(choices == input$enmchoice)]]@parameters$data, length(x@enms), algo.info, x@enms[[which(choices == input$enmchoice)]]@parameters$rep, PA)
      if(!is.null(x@enms[[which(choices == input$enmchoice)]]@parameters$sp.nb.origin)) {
        summary = rbind(summary, 
                        data.frame(Summary = x@enms[[which(choices == input$enmchoice)]]@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$enm.binary.info <- renderText({
      x@enms[[which(choices == input$enmchoice)]]@parameters$metric = switch(x@enms[[which(choices == input$enmchoice)]]@parameters$metric,
                                                                                             'Kappa' = 'maximizing Kappa',
                                                                                             'CCR' = 'maximizing correct proportion (CCR)',
                                                                                             'TSS' = 'maximizing sensitivity and specificity sum (TSS)',
                                                                                             'SES' = 'equalizing sensitivity and specificity',
                                                                                             'LW' = 'taking the minimum occurence prediction',
                                                                                             'ROC' = 'calculating the minimum ROC plot distance')
      text = paste('Binary map realized by', x@enms[[which(choices == input$enmchoice)]]@parameters$metric, 
                   'with a final threshold of',  round(x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold, digits = 3))
      text
    })
    output$enm.varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i])
        } else if (i == length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$enm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (x@enms[[which(choices == input$enmchoice)]]@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # ENM end #
  }
    
  shinyApp(ui, server)
  })
