#' @import methods
#' @importFrom sp Polygon Polygons SpatialPolygons bbox
#' @importFrom raster raster stack extract predict reclassify layerStats calc
#' @importFrom mgcv gam gam.control
#' @importFrom earth earth
#' @importFrom rpart rpart rpart.control
#' @importFrom gbm gbm
#' @importFrom randomForest randomForest
#' @importFrom dismo maxent
#' @importFrom nnet nnet
#' @importFrom e1071 svm
#' @importFrom grDevices heat.colors is.raster rainbow terrain.colors
#' @importFrom graphics arrows barplot legend
#' @importFrom stats aggregate.data.frame cor glm glm.control rbinom runif sd var
#' @importFrom utils lsf.str read.csv read.csv2 tail write.csv
NULL

setClass('SDM',
         representation(name = 'character',
                        projection = 'Raster',
                        binary = 'Raster',
                        evaluation = 'data.frame',
                        variable.importance = 'data.frame',
                        data = 'data.frame',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   projection = raster(),
                   binary = raster(),
                   evaluation = data.frame(),
                   variable.importance = data.frame(),
                   data = data.frame(),
                   parameters = data.frame()))

setMethod('print', 'SDM', function(x, ...) {
  cat('Object of class :', class(x)[1],'\n')
  cat('Name :', x@name, '\n')
  cat('Projections : ',names(x@projection),'\n')
  print(x@evaluation)
  print(x@variable.importance)
  if(inherits(x, 'Ensemble.SDM')) {
    cat('Uncertainty map :', names(x@uncertainty),'\n')
    print(x@algorithm.evaluation)
    print(x@algorithm.correlation)
  }
})
