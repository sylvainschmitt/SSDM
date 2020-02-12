#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

#' Project model into environment
#'
#'This is a collection of methods to project SDMs, ESDMs or SSDMs into the supplied environment. The function is used internally to calculate the input for the projection slot of .SDM classes but can also be used to project existing .SDM objects (see Details). 
#'
#' @param obj Object of class Algorithm.SDM, Ensemble.SDM or Stacked.SDM. Model(s) to be projected.
#' @param Env Raster stack. Updated environmental rasters to be used for projection.
#' @param ... Additional arguments for internal use.
#' @details  The function uses any S4 .SDM class object and a raster stack of environmental layers of the variables the model was trained with. 
#' @return Either returns the original .SDM object with updated projection slots or if minimal.outputs = TRUE only returns the projections as Raster* objects. Depending on the object class this may be: a raster (Algorithm.SDM), a raster stack (Ensemble.SDM), a biodiversity map/mean raster (Stacked.SDM).
#' @name project
#' @export
setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

#' @rdname project
#' @export
setMethod("project", "Algorithm.SDM", function(obj, Env, ...) {
  model = get_model(obj, ...)
  if(all(names(Env) %in% colnames(obj@data)[-c(1:3)])==FALSE){stop("Environmental layer names do not match the variables used for model training")}
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  if(length(factors)==0) factors <- NULL
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
    if(proj@data@max) proj = proj / proj@data@max
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
  if(all(names(Env) %in% colnames(obj@data)[-c(1:3)])==FALSE){stop("Environmental layer names do not match the variables used for model training")}
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
    if(proj@data@max) proj = proj / proj@data@max
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
  if(all(names(Env) %in% colnames(obj@data)[-c(1:3)])==FALSE){stop("Environmental layer names do not match the variables used for model training")}
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  if(length(factors)==0) factors <- NULL
  # project SDMs
  proj = suppressWarnings(lapply(models,FUN=function(x){raster::predict(Env, x, factors = factors)}))
  # rescaling
  proj = lapply(proj, FUN=function(x) reclassify(x, c(-Inf, 0, 0)))
  for(i in 1:length(models)){
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs should not be rescaled
      if(proj[[i]]@data@max>0) proj[[i]] = proj[[i]] / proj[[i]]@data@max
    names(proj[[i]]) = "Projection"
    obj@sdms[[i]]@projection = proj[[i]]
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs can't produce binary
      obj@sdms[[i]]@binary <- reclassify(proj[[i]], c(-Inf,obj@sdms[[i]]@evaluation$threshold,0, obj@sdms[[i]]@evaluation$threshold,Inf,1))
  }
  # sum SDMs (To do - enable minimal.outputs)
  #ensemble.args <- c(verbose=FALSE,ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")])
  sum.algo.ensemble <- do.call(ensemble, c(obj@sdms,list(ensemble.thresh=0,verbose=F,weight=obj@parameters[,which(names(obj@parameters)=="weight")])))
  obj@projection <- sum.algo.ensemble@projection
  obj@binary <- sum.algo.ensemble@binary
  obj@uncertainty <- sum.algo.ensemble@uncertainty
    
  return(obj)
})

#' @rdname project
#' @export
setMethod("project","Stacked.SDM",function(obj,Env,...){
  # get factors in Env
  if(all(names(Env) %in% colnames(obj@esdms[[1]]@data)[-c(1:3)])==FALSE){stop("Environmental layer names do not match the variables used for model training")}
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  if(length(factors)==0) factors <- NULL
  esdms <- list() # temporary object with ESDM by species (for use with get_model and predict)
  sum.algo.ensemble <- list() # temporary object for storing resulting ensembles
  species.names <- names(obj@esdms)
  for(j in 1:length(obj@esdms)){
    # get ESDMs by species
    esdms[[j]] <- lapply(obj@esdms[[j]]@sdms,FUN=get_model)
    # project SDMs
    proj = suppressWarnings(lapply(esdms[[j]],FUN=function(x){raster::predict(Env, x, factors = factors)}))
    # rescaling
    proj = lapply(proj, FUN=function(x) reclassify(x, c(-Inf, 0, 0)))
    for(i in 1:length(obj@esdms[[j]]@sdms)){
      if(all(obj@esdms[[j]]@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs should not be rescaled
        if(proj[[i]]@data@max>0) proj[[i]] = proj[[i]] / proj[[i]]@data@max # if zero is the maximum, then this leads to an NA map, which will propagate throughout the projections
      names(proj[[i]]) = "Projection"
      obj@esdms[[j]]@sdms[[i]]@projection = proj[[i]]
      if(all(obj@esdms[[j]]@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs can't produce binary
        obj@esdms[[j]]@sdms[[i]]@binary <- reclassify(proj[[i]], c(-Inf,obj@esdms[[j]]@sdms[[i]]@evaluation$threshold,0, obj@esdms[[j]]@sdms[[i]]@evaluation$threshold,Inf,1))
    }
    # sum SDMs (to do - use ensemble function with minimal.outputs)
    #ensemble.args <- list(verbose=FALSE,ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")])
    sum.algo.ensemble[[j]] <- do.call(ensemble, c(obj@esdms[[j]]@sdms,list(verbose=FALSE,ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")])))
    sum.algo.ensemble[[j]]@name <- species.names[j]
    obj@esdms[[j]]@projection <- sum.algo.ensemble[[j]]@projection
    obj@esdms[[j]]@binary <- sum.algo.ensemble[[j]]@binary
    obj@esdms[[j]]@uncertainty <- sum.algo.ensemble[[j]]@uncertainty
  } # end project ESDMs
  
  # stack ESDMs (To do - enable minimal.outputs)
  #stack.args <- list(verbose=FALSE,method=obj@parameters[,which(names(obj@parameters)=="method")])
  ensemble.stack <- do.call(stacking, c(sum.algo.ensemble,list(verbose=FALSE,method=obj@parameters[,which(names(obj@parameters)=="method")])))
  obj@diversity.map <- ensemble.stack@diversity.map
  obj@endemism.map <- ensemble.stack@endemism.map
  obj@uncertainty <- ensemble.stack@uncertainty
  return(obj)
})
