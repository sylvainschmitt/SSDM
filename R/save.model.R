#' @include Ensemble.SDM.R Stacked.SDM.R checkargs.R
#' @importFrom raster writeRaster
NULL

#' Save ensemble SDMs and SSDMs
#'
#' Allows to save S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' class objects.
#'
#' @param esdm Ensemble.SDM. Ensemble SDM to be saved.
#' @param stack Stacked.SDM. SSDM to be saved.
#' @param name character. Folder name of the model to save.
#' @param path character. Path to the directory chosen to save the SDM,
#'   by default the path to the current directory.
#' @param verbose logical. If set to true, allows the function to print text in the
#'   console.
#' @param GUI logical. Don't take that argument into account (parameter for the
#'   user interface).
#'
#' @return Nothing in R environment. Creates folders, tables and rasters
#'   associated to the SDM. Tables are in .csv and rasters in .grd/.gri.
#'
#' @seealso \code{\link{load.model}}
#'
#' @name save.model
#'
NULL

#' @rdname save.model
#' @export
setGeneric('save.esdm', function (esdm, name = strsplit(esdm@name, '.', fixed = TRUE)[[1]][1],
                                 path = getwd(), verbose = TRUE, GUI = FALSE) {return(standardGeneric('save.esdm'))})

#' @rdname save.model
#' @export
setMethod('save.esdm', 'Ensemble.SDM', function (esdm,
                                                        name = strsplit(esdm@name, '.Ensemble.SDM', fixed = TRUE)[[1]][1],
                                                        path = getwd(),
                                                        verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(esdm = esdm, name = name, path = path, verbose = verbose, GUI = GUI)

  if (verbose) {
    cat("Saving ensemble model results \n")
  }
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  dir.create(path = paste0(path, "/", name, "/Rasters"))
  dir.create(path = paste0(path, "/", name, "/Tables"))

  # Raster saving
  if (verbose) {
    cat("   rasters ...")
  }
  writeRaster(esdm@projection[[1]], paste0(path, "/", name, "/Rasters/Probability"),
              "GTiff", overwrite = TRUE)
  writeRaster(esdm@binary[[1]], paste0(path, "/", name, "/Rasters/Binary"),
              "GTiff", overwrite = TRUE)
  writeRaster(esdm@uncertainty, paste0(path, "/", name, "/Rasters/uncertainty"),
              "GTiff", overwrite = TRUE)
  if (verbose) {
    cat("saved \n")
  }

  # Tables saving
  if (verbose) {
    cat("   tables ...")
  }
  write.csv(esdm@evaluation, paste0(path, "/", name, "/Tables/esdmeval.csv"))
  write.csv(esdm@algorithm.evaluation, paste0(path, "/", name, "/Tables/AlgoEval.csv"))
  write.csv(esdm@algorithm.correlation, paste0(path, "/", name, "/Tables/AlgoCorr.csv"))
  write.csv(esdm@variable.importance, paste0(path, "/", name, "/Tables/VarImp.csv"))
  write.csv(esdm@data, paste0(path, "/", name, "/Tables/Data.csv"))
  write.csv(esdm@name, paste0(path, "/", name, "/Tables/Name.csv"))
  write.csv(esdm@parameters, paste0(path, "/", name, "/Tables/Parameters.csv"))
  if (verbose) {
    cat("saved \n \n")
  }
})

#' @rdname save.model
#' @export
setGeneric('save.stack', function (stack, name = 'Stack', path = getwd(), verbose = TRUE, GUI = FALSE) {return(standardGeneric('save.stack'))})

#' @rdname save.model
#' @export
setMethod('save.stack', 'Stacked.SDM', function (stack, name = 'Stack',
                                                                        path = getwd(),
                                                                        verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(stack = stack, name = name, path = path, verbose = verbose,
             GUI = GUI)

  if (verbose) {
    cat("Saving stack species model results \n")
  }
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  path = paste0(path, "/", name)
  dir.create(path = paste0(path, "/", "Stack"))
  dir.create(path = paste0(path, "/", "Species"))
  dir.create(path = paste0(path, "/", "Stack", "/Rasters"))
  dir.create(path = paste0(path, "/", "Stack", "/Tables"))

  # Raster saving
  if (verbose) {
    cat("   rasters ...")
  }
  writeRaster(stack@diversity.map, paste0(path, "/", "Stack", "/Rasters/Diversity"),
              "GTiff", overwrite = TRUE)
  writeRaster(stack@endemism.map, paste0(path, "/", "Stack", "/Rasters/Endemism"),
              "GTiff", overwrite = TRUE)
  writeRaster(stack@uncertainty, paste0(path, "/", "Stack", "/Rasters/uncertainty"),
              "GTiff", overwrite = TRUE)
  cat("saved \n")

  # Tables saving
  if (verbose) {
    cat("   tables ...")
  }
  write.csv(stack@evaluation, paste0(path, "/", "Stack", "/Tables/StackEval.csv"))
  write.csv(stack@algorithm.evaluation, paste0(path, "/", "Stack", "/Tables/AlgoEval.csv"))
  write.csv(stack@algorithm.correlation, paste0(path, "/", "Stack", "/Tables/AlgoCorr.csv"))
  write.csv(stack@variable.importance, paste0(path, "/", "Stack", "/Tables/VarImp.csv"))
  write.csv(stack@name, paste0(path, "/", "Stack", "/Tables/Name.csv"))
  write.csv(stack@parameters, paste0(path, "/", "Stack", "/Tables/Parameters.csv"))
  cat("saved \n\n")

  # ESDMs saving
  if (verbose) {
    cat("   esdms ... \n\n")
  }
  for (i in seq_len(length(stack@esdms))) {
    save.esdm(stack@esdms[[i]], path = paste0(path, "/", "Species"), verbose = verbose,
             GUI = GUI)
  }
  if (verbose) {
    cat("saved \n \n")
  }
})
