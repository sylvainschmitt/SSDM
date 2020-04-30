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
#' @param output.format character. If 'model' (default), the original .SDM object will be returned with updated projection slots. If 'rasters', the projected rasters will be returned as a list of rasters.
#' @param uncertainty logical. If set to TRUE, generates an uncertainty map. If output.format is 'model' an algorithm correlation matrix is additionally returned.
#' @param cores integer. Specify the number of CPU cores used to do the
#'  computing. You can use \code{\link[parallel]{detectCores}}) to automatically
#'  use all the available CPU cores.
#' @param minimal.memory logical. Only relevant if cores >1. If TRUE, only one model will be sent to each worker at a time, reducing used working memory.
#' @param tmp logical or character. If FALSE, no temporary rasters are written. If TRUE, temporary rasters are written to the „tmp“ directory of your R environment. If character, temporary rasters are written to a custom path. Very useful to reduce working memory consumption (use together with minimal.memory=TRUE for maximal effect).
#' But beware: Depending on number, resolution and extent of models, temporary files can take a lot of disk space. 
#' @param ... arguments for internal use (get_model), such as argument lists to be passed to the source functions (e.g. glm.args=list(test="AIC",singular.ok=FALSE)). See \code{\link[SSDM]{modelling}}, algorithm section for more details.
#' @details  The function uses any S4 .SDM class object and a raster stack of environmental layers of the variables the model was trained with. 
#' @return Either returns the original .SDM object with updated projection slots (default) or if output.format = 'rasters' only returns the projections as Raster* objects or a list thereof.
#' @name project
#' @export
setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

#' @rdname project
#' @export
setMethod("project", "Algorithm.SDM", function(obj, Env, output.format='model', ...) {
  if(!inherits(output.format,"character")||!output.format%in%c('model','rasters')){
    stop("output.format should be 'model' or 'rasters'.")
  }
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
  # proj = reclassify(proj, c(-Inf, 0, 0))
  # clamp projection into interval 0-1
  proj = reclassify(proj, c(-Inf, 0, 0,1,Inf,1))
  if(!all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    # normalization is bad
    # if(proj@data@max) proj = proj / proj@data@max
    names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
    obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                     obj@evaluation$threshold,Inf,1))
  
  if(output.format=="rasters"){
    return(list(projection=obj@projection,binary=obj@binary))
  } else if(output.format=="model"){
    return(obj)
  }

})

#' @rdname project
#' @export
setMethod("project", "MAXENT.SDM", function(obj, Env, output.format='model', ...) {
  model = get_model(obj, ...)
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
  # proj = reclassify(proj, c(-Inf, 0, 0))
  # clamp projection into interval 0-1
  proj = reclassify(proj, c(-Inf, 0, 0,1,Inf,1))
  if(!all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    # normalization is bad
    # if(proj@data@max) proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
    obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                     obj@evaluation$threshold,Inf,1))
  
  if(output.format=="rasters"){
    return(list(projection=obj@projection,binary=obj@binary))
  } else if(output.format=="model"){
    return(obj)
  }
})

#' @rdname project
#' @export
setMethod("project", "Ensemble.SDM", function(obj, Env, uncertainty=TRUE, output.format='model', SDM.projections=FALSE, cores=0, minimal.memory=FALSE,tmp=FALSE,...) {
  
  models <- obj@sdms
  if (cores > 0 && requireNamespace("parallel", quietly = TRUE)) {
    require(foreach)
    if ((parallel::detectCores()) < cores) {
      cores <- parallel::detectCores()
      warning(paste("It seems you attributed more cores than your CPU has! Automatic reduction to",cores, "cores."))
    }
    # minimal memory option (loop through model chunks) for using many cores at once
    if(minimal.memory){
      # create indices to split models into chunks
      chunks <- split(1:length(models),ceiling(seq_along(1:length(models))/cores))
      proj <- NULL
      for(k in 1:length(chunks)){
        model_chunk <- models[chunks[[k]]]
        
        cl <- parallel::makePSOCKcluster(length(model_chunk))
        doParallel::registerDoParallel(cl)
        proj_chunk <- foreach(model_chunk=itertools::isplitVector(model_chunk, chunks=length(cl)),.packages = c("raster","SSDM"),.verbose=FALSE) %dopar% {
          lapply(model_chunk,project,Env = Env)
        }
        
        proj_chunk <- unlist(proj_chunk)
    
    # save temporary rasters    
    if(!isFALSE(tmp)){
      if(isTRUE(tmp)){
        tmppath <- get("tmpdir", envir = .PkgEnv)
      }
      if(is.character(tmp)){
        tmppath <- tmp
        }
      if (!dir.exists(paste0(tmppath, "/.models"))){
        dir.create(paste0(tmppath, "/.models"))
      }
      for(i in 1:length(proj_chunk)){
        proj_chunk[[i]]@projection <- writeRaster(proj_chunk[[i]]@projection, filename=paste0(tmppath, "/.models/proba_",proj_chunk[[i]]@name,"-",chunks[[k]][i],gsub(" |:|-","", Sys.time())))
        proj_chunk[[i]]@binary <- writeRaster(proj_chunk[[i]]@binary, filename=paste0(tmppath, "/.models/bin_",proj_chunk[[i]]@name,"-",chunks[[k]][i],gsub(" |:|-","", Sys.time())))
      }
    } # tmp
      parallel::stopCluster(cl)
      proj <- c(proj,proj_chunk)
      rm(proj_chunk)
      gc(verbose=FALSE)
      } # k
    } else {
      # normal mode
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      proj <- foreach(models=itertools::isplitVector(models, chunks=cores),.packages = c("raster","SSDM","itertools"),.verbose=F) %dopar% lapply(models,FUN=function(x){project(x,Env)})
    
    parallel::stopCluster(cl)
    proj <- unlist(proj,recursive = FALSE)
    }
  } else {
    # sequential
    proj = lapply(models,FUN=function(x){project(x,Env)})
  }
  
  # ensemble SDMs
  sum.algo.ensemble <- do.call(ensemble, c(proj,list(ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")], verbose=F, SDM.projections=SDM.projections, uncertainty=uncertainty)))

  if(output.format=="model"){
    return(sum.algo.ensemble)  
  }  
  else if(output.format=="rasters"){
    projls <- list(projection=sum.algo.ensemble@projection, binary=sum.algo.ensemble@binary)
    if(uncertainty){projls <- c(projls,uncertainty=sum.algo.ensemble@uncertainty)}
    if(SDM.projections){projls <- c(projls, list(sdms=lapply(sum.algo.ensemble@sdms, function(x) list(projection=x@projection,binary=x@binary))))}
    return(projls)
  }
  # clean up tmp -- does not make sense if individual projections should be kept
  # if(isTRUE(tmp)|is.character(tmp)){
  #   unlink(paste0(tmppath,"/.models"), recursive = TRUE, force = TRUE)
  # }
})

#' @rdname project
#' @export
setMethod("project","Stacked.SDM",function(obj,Env,method=NULL, uncertainty=TRUE, output.format='model', SDM.projections=FALSE, cores=0, minimal.memory=FALSE, tmp=FALSE, ...){
  
  esdms <- obj@esdms
  
  esdms_proj <- lapply(esdms,function(x){
    project(obj=x,Env=Env,uncertainty=uncertainty,cores=cores,SDM.projections=SDM.projections,minimal.memory=minimal.memory,tmp=tmp)
  })
  # transfer model names, so they are not used as argument names
   for(i in 1:length(esdms_proj)){
    esdms_proj[[i]]@name <- names(esdms_proj)[i]
    names(esdms_proj)[i] <- ""
   }
  # for(j in 1:length(obj@esdms)){
  #   
  #   # ensemble SDMs
  #   sum.algo.ensemble[[j]] <- do.call(ensemble, c(obj@esdms[[j]]@sdms,list(ensemble.thresh=0,weight=obj@parameters[,which(names(obj@parameters)=="weight")], verbose=F, SDM.projections=SDM.projections, uncertainty=uncertainty)))
  #   sum.algo.ensemble[[j]]@name <- species.names[j]
  #   # obj@esdms[[j]]@projection <- sum.algo.ensemble[[j]]@projection
  #   # obj@esdms[[j]]@binary <- sum.algo.ensemble[[j]]@binary
  #   # obj@esdms[[j]]@uncertainty <- sum.algo.ensemble[[j]]@uncertainty
  # } # end project ESDMs
  # 
  # stack ESDMs 
  if(is.null(method)){method <- obj@parameters[,which(names(obj@parameters)=="method")]}
  ensemble.stack <- do.call(stacking, c(esdms_proj,list(verbose=FALSE,method=method, uncertainty=uncertainty)))
  # obj@diversity.map <- ensemble.stack@diversity.map
  # obj@endemism.map <- ensemble.stack@endemism.map
  # obj@uncertainty <- ensemble.stack@uncertainty
  if(output.format=="model"){
    return(ensemble.stack)
  } else if(output.format=="rasters"){
    projls <- list(diversity.map=ensemble.stack@diversity.map, endemism.map=ensemble.stack@endemism.map, esdms=lapply(ensemble.stack@esdms, function(x) list(projection=x@projection,binary=x@binary)))
    if(uncertainty){
      projls <- c(projls,uncertainty=ensemble.stack@uncertainty)
      projls$esdms <- c(projls$esdms,list(lapply(ensemble.stack@esdms, function(x) list(uncertainty=x@uncertainty))))
    }
    if(SDM.projections){projls$esdms <- c(projls$esdms, list(sdms=lapply(ensemble.stack@esdms, function(x) lapply(x@sdms, function(y) list(projection=y@projection,binary=y@binary)))))}
    return(projls)
  }
  
})
