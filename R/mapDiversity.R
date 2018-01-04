#' @include Ensemble.SDM.R checkargs.R
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPoints bbox
#' @importFrom raster raster stack reclassify mask calc overlay values rasterize rasterToPoints values<-
#' @importFrom stats lm
NULL

#' Map Diversity
#'
#' Methods for Stacked.SDM or SSDM to map diversity and communities composition.
#'
#'@param obj Stacked.SDM. SSDM to map diversity with.
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random Bernoulli (\strong{Bernoulli}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'@param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}). Needed only for stacking
#'  method using probability ranking from richness (\strong{PRR}).
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'@param ... other arguments pass to the method.
#'
#'@return a list with a diversity map and eventually ESDMs for stacking method
#'  using probability ranking from richness (\strong{PPR}).
#'
#'@details \strong{Methods:} Choice of the method used to compute the local
#'  species richness map (see Calabrez et al. (2014) and D'Amen et al (2015) for
#'  more informations, see reference below): \describe{\item{pSSDM}{sum
#'  probabilities of habitat suitability maps}\item{Bernoulli}{draw repeatedly
#'  from a Bernoulli distribution}\item{bSSDM}{sum the binary map obtained with
#'  the thresholding (depending on the metric of the
#'  ESDM).}\item{MaximumLikelihood}{adjust species richness of the model by
#'  linear regression}\item{PRR.MEM}{model richness with a macroecological model
#'  (MEM) and adjust each ESDM binary map by ranking habitat suitability and
#'  keeping as much as predicted richness of the MEM}\item{PRR.pSSDM}{model
#'  richness with a pSSDM and adjust each ESDM binary map by ranking habitat
#'  suitability and keeping as much as predicted richness of the pSSDM}}
#'
#'@examples
#'
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Diversity mapping
#' mapDiversity(SSDM, mathod = 'pSSDM')
#'
#' }
#'
#'@seealso \code{\link{stacking}} to build SSDMs.
#'
#'@references M. D'Amen, A. Dubuis, R. F. Fernandes, J. Pottier, L. Pelissier, &
#'  A Guisan (2015) "Using species richness and functional traits prediction to
#'  constrain assemblage predicitions from stacked species distribution models"
#'  \emph{Journal of Biogeography} 42(7):1255-1266
#'  \url{http://doc.rero.ch/record/235561/files/pel_usr.pdf}
#'
#'  J.M. Calabrese, G. Certain, C.  Kraan, & C.F. Dormann (2014) "Stacking
#'  species distribution  models  and  adjusting  bias  by linking them to
#'  macroecological models." \emph{Global Ecology and Biogeography} 23:99-112
#'  \url{http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/calabrese2013globalecolbiogeogr.pdf}
#'
#' @name mapDiversity
NULL

#' @rdname mapDiversity
#' @export
setGeneric('mapDiversity', function(obj, ...) {return(standardGeneric('mapDiversity'))})

#' @rdname mapDiversity
#' @export
setMethod("mapDiversity", "Stacked.SDM", function(obj, method, rep.B = 1000,
                                                  verbose = TRUE, Env = NULL,
                                                  ...){
  # Check arguments
  .checkargs(stack = obj, method = method, rep.B = rep.B, verbose = verbose)

  enms <- NULL # Preparing enms slot for PPR methods
  diversity.map <- reclassify(obj@enms[[1]]@projection[[1]], c(-Inf,Inf, 0))

  # Useless datacheck to prevent bugs to remove after debugging
  for (i in seq_len(length(obj@enms))) {
    if (!inherits(obj@enms[[i]]@projection, "RasterLayer")) {
      if (verbose) {
        cat("Error", obj@enms[[i]]@name, "is not a raster but a",
            class(obj@enms[[i]]@projection)[1],
            ".\nIt will be removed for the stacking")
      }
      obj@enms[[i]] <- NULL
    }
  }

  if (method == "bSSDM") {
    # Threshold and sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by thresholding and then summing. \n")
    diversity.map <- sum(stack(lapply(obj@enms, function(x)
      reclassify(x@projection,
                 c(-Inf,x@evaluation$threshold,0,x@evaluation$threshold,Inf,1))
      )))
  }

  if (method == "pSSDM") {
    # Individual probabilities sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by summing individual probabilities. \n")
    diversity.map <- sum(stack(lapply(obj@enms, function(x) x@projection)))
  }

  if (method == "Bernoulli") {
    # Random Bernoulli distribution (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by drawing repeatedly from a Bernoulli distribution. \n")
    proba <- stack(lapply(obj@enms, function(x) x@projection))
    diversity.map <- calc(proba, fun = function(...) {
      x <- c(...)
      x[is.na(x)] <- 0
      return(rbinom(lengths(x), rep.B, x))
    }, forcefun = TRUE)
    diversity.map <- sum(diversity.map)/length(enms)/rep.B
  }

  if (method == "MaximumLikelihood") {
    # Maximum likelihood (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by maximum likelihood adjustment. \n")
    diversity.map <- mapDiversity(obj, method = 'bSSDM',
                                  verbose = FALSE)$diversity.map
    Richness <- .richness(obj)
    SSDM_Richness <- values(mask(diversity.map, Richness))
    SSDM_Richness <- SSDM_Richness[-which(is.na(SSDM_Richness))]
    Richness <- values(Richness)
    Richness <- Richness[-which(is.na(Richness))]
    fit <- lm(Richness ~ SSDM_Richness)
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    diversity.map <- a + b * diversity.map
  }

  if (method == "PRR.MEM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from MEM. \n")
    diversity.map <- .MEM(obj,Env)@projection
    enms <- .PRR(obj, diversity.map)
  }

  if (method == "PRR.pSSDM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from pSSDM. \n")
    diversity.map <- mapDiversity(obj, method = 'pSSDM',
                                  verbose = FALSE)$diversity.map
    enms <- .PRR(obj, diversity.map)
  }

  return(list(
    diversity.map = diversity.map,
    enms = enms
  ))
})

##### Internals ####

.richness <- function(obj){
  Richness <- reclassify(obj@enms[[1]]@projection, c(-Inf, Inf, 0))
  for (i in seq_len(length(obj@enms)))
    Richness <- Richness + rasterize(
      SpatialPoints(obj@enms[[i]]@data[1:2]),
      Richness, field = obj@enms[[i]]@data$Presence,
      background = 0)
  if (all(values(Richness) %in% c(0, 1, NA)))
    stop("Observed Richness is always equal to 1, modelled richness can't be adjusted !")
  return(Richness)
}

.MEM <- function(obj, Env){
  occ <- data.frame(rasterToPoints(.richness(obj), function(x) x > 0))
  ensemble_modelling(algorithms = unlist(
    strsplit(obj@enms[[1]]@parameters$algorithms,
             ".", fixed = TRUE))[-1],
    Occurrences = occ, Env = Env, Xcol = "x",
    Ycol = "y", Pcol = "layer", rep = obj@enms[[1]]@parameters$rep,
    name = "MEM", cv = obj@enms[[1]]@parameters$cv,
    cv.param = as.numeric(unlist(
      strsplit(obj@enms[[1]]@parameters$cv.param,
               "|", fixed = TRUE))[-1]),
    metric = obj@enms[[1]]@parameters$metric,
    axes.metric = obj@enms[[1]]@parameters$axes.metric,
    ensemble.metric = unlist(
      strsplit(obj@enms[[1]]@parameters$ensemble.metric,
               ".", fixed = TRUE))[-1],
    ensemble.thresh = as.numeric(unlist(
      strsplit(obj@enms[[1]]@parameters$ensemble.thresh,
               "|", fixed = TRUE))[-1]),
    uncertainty = FALSE,
    weight = as.logical(obj@enms[[1]]@parameters$weight),
    verbose = FALSE)
}

.PRR <- function(obj, Richness){
  # Readjust each enm binary map
  richnesses <- values(Richness)
  names(richnesses) <- seq_len(length(richnesses))
  richnesses <- as.list(richnesses)
  probabilities <- lapply(lapply(obj@enms, FUN = slot, name = "projection"),
                          values)
  probabilities <- lapply(probabilities, function(x) {
    names(x) <- rep(seq_len(length(probabilities[[1]])))
    return(x)
  })
  probabilities <- lapply(probabilities, `[`, names(probabilities[[1]]))
  probabilities <- apply(do.call(rbind, probabilities), 2, as.list)
  binaries <- lapply(lapply(obj@enms, FUN = slot, name = "binary"),
                     values)
  binaries <- lapply(binaries, function(x) {
    names(x) <- rep(seq_len(length(binaries[[1]])))
    return(x)
  })
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- mapply(function(rich, probability, binary) {
    if (!is.na(rich)) {
      ord <- order(unlist(probability), decreasing = TRUE)
      binary <- unlist(binary)
      if (length(ord) <= rich) {
        binary[ord] <- 1
      } else {
        binary[ord[1:rich]] <- 1
        binary[ord[rich + seq_len(length(ord))]] <- 0
      }
      binary <- as.list(binary)
    }
    return(binary)
  }, rich = richnesses, probability = probabilities, binary = binaries,
  SIMPLIFY = FALSE)
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- lapply(binaries, unlist)
  binaries <- lapply(binaries, unname)
  mapply(function(enm, binary) {
    values(enm@binary) <- binary
    return(enm)
  }, enm = obj@enms, binary = binaries, SIMPLIFY = FALSE)
}
