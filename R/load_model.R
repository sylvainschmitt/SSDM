#' @include Ensemble.SDM.R Stacked.SDM.R
#' @importFrom shiny incProgress
#' @importFrom raster raster stack reclassify
NULL

#' Load ensemble SDMs and SSDMs.
#'
#' Load S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' objects saved with their respective save function.
#'
#' @param name character. Name of the folder containing the model to be
#'   loaded.
#' @param path character. Path to the directory containing the model to be
#'   loaded, by default the path to the current directory.
#' @param GUI logical. Do not take this argument into account (parameter for the
#'   user interface).
#'
#' @return The corresponding SDM object.
#'
#' @seealso \code{\link{save.model}}
#'
#' @name load.model
#'
NULL

#' @rdname load.model
#' @export
load_esdm <- function(name, path = getwd()) {
  path <- paste0(path, "/", name)
  a <- try(read.csv(paste0(path, "/Tables/AlgoCorr.csv"), row.names = 1))
  if (inherits(a, "try-error")) {
    cat("Algorithm correlation table empty ! \n")
    a <- data.frame()
  }
  b <- try(raster(paste0(path, "/Rasters/Binary.tif")))
  if (inherits(b, "try-error")) {
    projection <- raster(paste0(path, "/Rasters/Probability.tif"))
    evaluation <- read.csv(paste0(path, "/Tables/ESDMeval.csv"), row.names = 1)
    b <- reclassify(projection, c(-Inf, evaluation$threshold, 0, evaluation$threshold,
                                  Inf, 1))
  }
  esdm <- Ensemble.SDM(name = as.character(read.csv(paste0(path, "/Tables/Name.csv"))[1,2]),
                      projection = raster(paste0(path, "/Rasters/Probability.tif")),
                      binary = b,
                      uncertainty = try(raster(paste0(path, "/Rasters/uncertainty.tif"))),
                      evaluation = read.csv(paste0(path, "/Tables/ESDMeval.csv"), row.names = 1),
                      algorithm.evaluation = read.csv(paste0(path, "/Tables/AlgoEval.csv"),
                                                      row.names = 1),
                      algorithm.correlation = a,
                      data = read.csv(paste0(path, "/Tables/Data.csv"), row.names = 1),
                      variable.importance = read.csv(paste0(path, "/Tables/VarImp.csv"), row.names = 1),
                      parameters = read.csv(paste0(path, "/Tables/Parameters.csv"), row.names = 1, colClasses = "character"))
  return(esdm)
}

#' @rdname load.model
#' @export
load_stack <- function(name = "Stack", path = getwd(), GUI = FALSE) {
  path <- paste0(path, "/", name)
  a <- try(read.csv(paste0(path, "/Stack/Tables/AlgoCorr.csv"), row.names = 1))
  if (inherits(a, "try-error")) {
    cat("Algorithm correlation table empty ! \n")
    a <- data.frame()
  }
  stack <- Stacked.SDM(name = as.character(read.csv(paste0(path, "/Stack/Tables/Name.csv"))[1,
                                                                                            2]), diversity.map = raster(paste0(path, "/Stack/Rasters/Diversity.tif")),
                       endemism.map = raster(paste0(path, "/Stack/Rasters/Endemism.tif")),
                       uncertainty = raster(paste0(path, "/Stack/Rasters/uncertainty.tif")),
                       evaluation = read.csv(paste0(path, "/Stack/Tables/StackEval.csv"),
                                             row.names = 1), variable.importance = read.csv(paste0(path,
                                                                                                   "/Stack/Tables/VarImp.csv"), row.names = 1), algorithm.correlation = a,
                       algorithm.evaluation = read.csv(paste0(path, "/Stack/Tables/AlgoEval.csv"),
                                                       row.names = 1), esdms = list(), parameters = read.csv(paste0(path,
                                                                                                                   "/Stack/Tables/Parameters.csv"), row.names = 1, colClasses = "character"))
  esdms <- list.dirs(paste0(path, "/Species"), recursive = FALSE, full.names = FALSE)
  if (GUI) {
    incProgress(1/(length(esdms) + 1), detail = "stack main results")
  }
  for (i in seq_len(length(esdms))) {
    cat(esdms[i])
    esdm <- load_esdm(esdms[i], path = paste0(path, "/Species"))
    stack@esdms[[esdm@name]] <- esdm
    if (GUI) {
      incProgress(1/(length(esdms) + 1), detail = esdm@name)
    }
  }
  return(stack)
}
