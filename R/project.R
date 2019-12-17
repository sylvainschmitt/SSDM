#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom itertools isplitVector
NULL

#' Project model into environment
#'
#'This is a collection of methods to project SDMs, ESDMs or SSDMs into the supplied environment. The function is used internally to calculate the input for the projection slot of .SDM classes but can also be used to project existing .SDM objects (see Details). 
#'
#' @param obj Object of class Algorithm.SDM, Ensemble.SDM or Stacked.SDM. Model(s) to be projected.
#' @param Env Raster stack. Updated environmental rasters to be used for projection.
#' @param method character. Define the method used to create the local species
#'  richness map (for details see \code{\link[SSDM]{stack_modelling}}). If NULL (default), the method used for building the SSDM is used.
#' @param SDM.projections logical. If FALSE (default), the projections of the Algorithm.SDMs will not be returned (only applies to Ensemble.SDMs and Stack.SDMs).
#' @param update.projections logical. If TRUE (default), the original .SDM object will be returned with the difference, that the projected rasters will replace the existing ones inside the respective projection slots. If FALSE, the projected rasters will be returned as a separate raster object (for single models) or a list of rasters (for ensembles and stacks).
#' @param uncertainty logical. If set to TRUE, generates an uncertainty map and if update.projection is TRUE
#'  additionally an algorithm correlation matrix.
#' @param cores integer. Specify the number of CPU cores used to do the
#'  computing. You can use \code{\link[parallel]{detectCores}}) to automatically
#'  use all the available CPU cores.
#' @param ... Additional arguments for internal use.
#' @details  The function uses any S4 .SDM class object and a raster stack of environmental layers of the variables the model was trained with. 
#' @return Either returns the original .SDM object with updated projection slots (default) or if update.projections = FALSE only returns the projections as Raster* objects or a list thereof.
#' @name project
#' @export
setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

#' @rdname project
#' @export
setMethod("project", "Algorithm.SDM", function(obj, Env, update.projections=TRUE,...) {
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
  
  if(!update.projections){
    return(obj@projection)
  } else {
    return(obj)
  }

})

#' @rdname project
#' @export
setMethod("project", "MAXENT.SDM", function(obj, Env, update.projections=TRUE, ...) {
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
  if(!update.projections){
    return(obj@projection)
  } else{
    return(obj)
  }
})

#' @rdname project
#' @export
setMethod("project", "Ensemble.SDM", function(obj, Env, uncertainty=TRUE, update.projections=TRUE, SDM.projections=FALSE, cores=0, ...) {
  models = lapply(obj@sdms,FUN=get_model)
  if(all(names(Env) %in% colnames(obj@data)[-c(1:3)])==FALSE){stop("Environmental layer names do not match the variables used for model training")}
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  if(length(factors)==0) factors <- NULL
  # project SDMs
  
  if (cores > 0 && requireNamespace("parallel", quietly = TRUE)) {
    if ((parallel::detectCores() - 1) < cores) {
      cores <- parallel::detectCores()-1
      warning(paste("It seems you attributed more cores than your CPU has! Automatic reduction to",cores, "cores."))
    }
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    proj <- foreach::foreach(models=itertools::isplitVector(models, chunks=cores),.packages = c("raster","itertools"),.verbose=F) %dopar% lapply(models,FUN=function(x){
      p = suppressWarnings(predict(object=Env,model=x,factors=factors))
      # rescale
      p = reclassify(p, c(-Inf, 0, 0))
      return(p)
    })
    proj <- unlist(proj,recursive = FALSE)
    parallel::stopCluster(cl)
  } else {
    proj = suppressWarnings(lapply(models,FUN=function(x){raster::predict(Env, x, factors = factors)}))
    # rescaling
    proj = lapply(proj, FUN=function(x) reclassify(x, c(-Inf, 0, 0)))
  }
  for(i in 1:length(models)){
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs should not be rescaled
      if(proj[[i]]@data@max>0) proj[[i]] = proj[[i]] / proj[[i]]@data@max
    names(proj[[i]]) = paste("Projection",obj@sdms[[i]]@name)
    obj@sdms[[i]]@projection = proj[[i]]
    if(all(obj@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs can't produce binary
      obj@sdms[[i]]@binary <- reclassify(proj[[i]], c(-Inf,obj@sdms[[i]]@evaluation$threshold,0, obj@sdms[[i]]@evaluation$threshold,Inf,1))
  }
  # ensemble SDMs
  sum.algo.ensemble <- do.call(ensemble, c(obj@sdms,list(ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")], verbose=F, SDM.projections=SDM.projections, uncertainty=uncertainty)))

  if(update.projections){
    return(sum.algo.ensemble)  
  }  
  else {
    projls <- list(projection=sum.algo.ensemble@projection, binary=sum.algo.ensemble@binary)
    if(uncertainty){projls <- c(projls,uncertainty=sum.algo.ensemble@uncertainty)}
    if(SDM.projections){projls <- c(projls, list(sdms=lapply(sum.algo.ensemble@sdms, function(x) list(projection=x@projection,binary=x@binary))))}
    return(projls)
  }
})

#' @rdname project
#' @export
setMethod("project","Stacked.SDM",function(obj,Env,method=NULL, uncertainty=TRUE, update.projections=TRUE, SDM.projections=FALSE, cores=0, ...){
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
  
  if ((parallel::detectCores() - 1) < cores) {
    cores <- parallel::detectCores()-1
    warning(paste("It seems you attributed more cores than your CPU has! Automatic reduction to",cores, "cores."))
  }
  
  for(j in 1:length(obj@esdms)){
    # get ESDMs by species
    esdms[[j]] <- lapply(obj@esdms[[j]]@sdms,FUN=get_model)
    
    if (cores > 0 && requireNamespace("parallel", quietly = TRUE)) {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      # project SDMs
      models <- esdms[[j]]
      proj <- foreach::foreach(models=itertools::isplitVector(models, chunks=cores),.packages = c("raster","itertools"),.verbose=F) %dopar% lapply(models,FUN=function(x){
        p = suppressWarnings(predict(object=Env,model=x,factors=factors))
        # rescale
        p = reclassify(p, c(-Inf, 0, 0))
        return(p)
      })
      proj <- unlist(proj,recursive = FALSE)
      parallel::stopCluster(cl)
      } else {
        # project SDMs
        proj = suppressWarnings(lapply(esdms[[j]],FUN=function(x){raster::predict(Env, x, factors = factors)}))
        # rescaling
        proj = lapply(proj, FUN=function(x) reclassify(x, c(-Inf, 0, 0)))
      }
    for(i in 1:length(obj@esdms[[j]]@sdms)){
      if(all(obj@esdms[[j]]@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs should not be rescaled
        if(proj[[i]]@data@max>0) proj[[i]] = proj[[i]] / proj[[i]]@data@max # if zero is the maximum, then this leads to an NA map, which will propagate throughout the projections
      names(proj[[i]]) = "Projection"
      obj@esdms[[j]]@sdms[[i]]@projection = proj[[i]]
      if(all(obj@esdms[[j]]@sdms[[i]]@data$Presence %in% c(0,1))) # MEMs can't produce binary
        obj@esdms[[j]]@sdms[[i]]@binary <- reclassify(proj[[i]], c(-Inf,obj@esdms[[j]]@sdms[[i]]@evaluation$threshold,0, obj@esdms[[j]]@sdms[[i]]@evaluation$threshold,Inf,1))
    }
    # ensemble SDMs
    sum.algo.ensemble[[j]] <- do.call(ensemble, c(obj@esdms[[j]]@sdms,list(ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")], verbose=F, SDM.projections=SDM.projections, uncertainty=uncertainty)))
    sum.algo.ensemble[[j]]@name <- species.names[j]
    # obj@esdms[[j]]@projection <- sum.algo.ensemble[[j]]@projection
    # obj@esdms[[j]]@binary <- sum.algo.ensemble[[j]]@binary
    # obj@esdms[[j]]@uncertainty <- sum.algo.ensemble[[j]]@uncertainty
  } # end project ESDMs
  
  # stack ESDMs 
  if(is.null(method)){method <- obj@parameters[,which(names(obj@parameters)=="method")]}
  ensemble.stack <- do.call(stacking, c(sum.algo.ensemble,list(verbose=FALSE,method=method, uncertainty=uncertainty)))
  # obj@diversity.map <- ensemble.stack@diversity.map
  # obj@endemism.map <- ensemble.stack@endemism.map
  # obj@uncertainty <- ensemble.stack@uncertainty
  if(update.projections){
    return(ensemble.stack)
  }
  else {
    projls <- list(diversity.map=ensemble.stack@diversity.map, endemism.map=ensemble.stack@endemism.map, esdms=lapply(ensemble.stack@esdms, function(x) list(projection=x@projection,binary=x@binary)))
    if(uncertainty){
      projls <- c(projls,uncertainty=ensemble.stack@uncertainty)
      projls$esdms <- c(projls$esdms,list(lapply(ensemble.stack@esdms, function(x) list(uncertainty=x@uncertainty))))
    }
    if(SDM.projections){projls$esdms <- c(projls$esdms, list(sdms=lapply(ensemble.stack@esdms, function(x) lapply(x@sdms, function(y) list(projection=y@projection,binary=y@binary)))))}
    return(projls)
  }
  
})
