#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom sp Polygon Polygons SpatialPolygons bbox
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

setGeneric('PA.select', function(obj, Env, ...) {return(standardGeneric('PA.select'))})

setMethod('PA.select', "Algorithm.SDM", function(obj, Env, PA = NULL, verbose = TRUE) {
  if (is.null(PA)) {
    PA = get_PA(obj)
    obj@parameters$PA = 'default'
  } else {
    obj@parameters$PA = paste0(as.character(PA$nb),'.',as.character(PA$strat))
  }

  # Mask defining
  if (PA$strat == 'disk') {
    if(verbose) {cat('   disk selection \n')}
    circles = list()
    for (i in seq_len(length(obj@data$X))) {
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
    if(verbose) {cat('   random selection \n')}
    Mask = Env[[1]]
  }

  # Pseudo-Absences selection
  data.PA = data.frame(matrix(nrow = 0, ncol = 2))
  names(data.PA) = c('X','Y')
  if(PA$nb < 100) {nb = PA$nb*PA$nb} else {nb = 1000}
  while (length(data.PA[,1]) < PA$nb) {
    X = runif(nb, min = bbox(Mask)[1,1], max = bbox(Mask)[1,2])
    Y = runif(nb, min = bbox(Mask)[2,1],max = bbox(Mask)[2,2])
    points = data.frame(X = X, Y = Y)
    check = extract(Mask, points)
    if(length(is.na(check)) > 0)
      points = points[-which(is.na(check)),]
    data.PA = rbind(data.PA, points)
  }
  data.PA = data.PA[1:PA$nb,]
  data.PA$Presence = 0
  # A bug can appear here if obj@data coordinates are consider as factors,
  # it can happen if occurences are not well loaded and decimal is not well defined
  obj@data = rbind(obj@data, data.PA)

  return(obj)})
