#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

#' Project model into a new environment
#'
#'This is a method to project existing SDMs, ESDMs or SSDMs to a new environment. The function uses any S4 .SDM class object and either returns the object with updated projection slots or only returns the projection as rasters (if minimal.outputs = TRUE)
#'
#' @param obj Object of class Algorithm.SDM, Ensemble.SDM or Stacked.SDM. Model(s) to be projected.
#' @param Env Raster stack. Updated environmental rasters to be used for projection.
#' @param ... 
#'
#' @return A raster (Algorithm.SDM), raster stack (Ensemble.SDM), biodiversity map/mean raster (Stacked.SDM)
#' @name project
#' @export
setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

#' @rdname project
#' @export
setMethod("project", "Algorithm.SDM", function(obj, Env, ...) {
  model = get_model(obj, ...)
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  proj = suppressWarnings(raster::predict(Env, model, factors = factors))
  # proj = suppressWarnings(raster::predict(Env, model, fun = function(model,
  #                                                                    x) {
  #   x = as.data.frame(x)
  #   for (i in seq_len(length(Env@layers))) {
  #     if (Env[[i]]@data@isfactor) {
  #       x[, i] = as.factor(x[, i])
  #       x[, i] = droplevels(x[, i])
  #       levels(x[, i]) = Env[[i]]@data@attributes[[1]]$ID
  #     }
  #   }
  #   return(predict(model, x))
  # }))
  # Rescaling projection
  proj = reclassify(proj, c(-Inf, 0, 0))
  if(all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
  obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                   obj@evaluation$threshold,Inf,1))
  return(obj)
})

#' @rdname project
#' @export
setMethod("project", "MAXENT.SDM", function(obj, Env, ...) {
  model = get_model(obj, Env, ...)
  proj = raster::predict(Env, model, fun = function(model, x) {
    x = as.data.frame(x)
    for (i in seq_len(length(Env@layers))) {
      if (Env[[i]]@data@isfactor) {
        x[, i] = as.factor(x[, i])
        x[, i] = droplevels(x[, i])
        levels(x[, i]) = Env[[i]]@data@attributes[[1]]$ID
      }
    }
    return(predict(model, x))
  })
  # Rescaling projection
  proj = reclassify(proj, c(-Inf, 0, 0))
  if(!all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
    obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                     obj@evaluation$threshold,Inf,1))
  return(obj)
})

#' @rdname project
#' @export
setMethod("project", "Ensemble.SDM", function(obj, Env, ...) {
  models = lapply(obj@sdms,FUN=get_model)
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  # project SDMs
  proj = suppressWarnings(lapply(models,FUN=function(x){raster::predict(Env, x, factors = factors)}))
  # rescaling
  proj = lapply(proj, FUN=function(x) reclassify(x, c(-Inf, 0, 0)))
  for(i in 1:length(models)){
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs should not be rescaled
      proj[[i]] = proj[[i]] / proj[[i]]@data@max
    names(proj[[i]]) = "Projection"
    obj@sdms[[i]]@projection = proj[[i]]
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs can't produce binary
      obj@sdms[[i]]@binary <- reclassify(proj[[i]], c(-Inf,obj@sdms[[i]]@evaluation$threshold,0, obj@evaluation$threshold,Inf,1))
  }
  # sum SDMs (use ensemble function with minimal.outputs)
  ensemble.args <- list(verbose=FALSE)
  sum.algo.ensemble <- do.call(ensemble, c(obj@sdms,ensemble.args))
  obj@projection <- sum.algo.ensemble@projection
    
  return(obj)
})
