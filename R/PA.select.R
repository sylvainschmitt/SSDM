#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom sf st_as_sf st_buffer
#' @importFrom raster raster stack extract predict reclassify layerStats calc Which xyFromCell crs
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
    circles <- st_as_sf(obj@data, coords = c("X", "Y"))
    circles <- st_buffer(circles, 2/60)
    Mask = mask(Env[[1]], circles)
  }
  if(PA$strat=='geobuffer'){
    if(!is.null(PA$dist)){
      dist_m <- PA$dist
    } else {dist_m <- 1000}
    if (!requireNamespace("geobuffer", quietly = TRUE))
      stop("Package \"geobuffer\" needed for this function to work. Please install it.",
           call. = FALSE)
    if (!requireNamespace("fasterize", quietly = TRUE))
      stop("Package \"fasterize\" needed for this function to work. Please install it.",
           call. = FALSE)
    buf_pts <- geobuffer::geobuffer_pts(xy = cbind(obj@data$X,obj@data$Y),
                                        dist_m = dist_m,crs=as.character(crs(Env[[1]])),
                                        output = "sf")
    buf_ras <- fasterize::fasterize(sf = buf_pts,
                                    raster = Env[[1]])
    Mask <- mask(Env[[1]], buf_ras, inverse = TRUE)
  } else {
    if(verbose) {cat('   random selection \n')}
    Mask = Env[[1]]
  }


  # Pseudo-Absences selection
  mask_inds <- Which(!is.na(Mask),cells=TRUE)
  if(PA$nb > length(mask_inds)){
    PA$nb <- length(mask_inds)
    if(verbose) {cat('Number of requested pseudo-absences exceeds number of cells, lowered to number of cells \n')}
  }

  abs_inds <- sample(mask_inds, PA$nb)
  abs <- xyFromCell(Mask,abs_inds)
  data.PA <- data.frame(X=abs[,1],Y=abs[,2],Presence=0)

  # old code follows
  # data.PA = data.frame(matrix(nrow = 0, ncol = 2))
  # names(data.PA) = c('X','Y')
  # #if(PA$nb < 100) {nb = PA$nb*PA$nb} else {nb = 1000}
  # while (length(data.PA[,1]) < PA$nb) {
  #   X = runif(nb, min = bbox(Mask)[1,1], max = bbox(Mask)[1,2])
  #   Y = runif(nb, min = bbox(Mask)[2,1],max = bbox(Mask)[2,2])
  #   points = data.frame(X = X, Y = Y)
  #   check = extract(Mask, points)
  #   if(length(is.na(check)) > 0)
  #     points = points[-which(is.na(check)),]
  #   data.PA = rbind(data.PA, points)
  # }
  # data.PA = data.PA[1:PA$nb,]
  # data.PA$Presence = 0
  # A bug can appear here if obj@data coordinates are consider as factors,
  # it can happen if occurences are not well loaded and decimal is not well defined
  obj@data = rbind(obj@data, data.PA)

  return(obj)})
