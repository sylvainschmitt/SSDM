#' @include Ensemble.SDM.R checkargs.R
#' @importFrom raster raster stack reclassify mask calc overlay values rasterize rasterToPoints values<- Which setValues
#' @importFrom stats lm optim
#' @importFrom poibin dpoibin
#' @importFrom sf st_as_sf
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
#'  species richness map (see Calabrese et al. (2014) and D'Amen et al (2015) for
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
#'  \url{https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12102}
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

  esdms <- NULL # Preparing esdms slot for PPR methods
  diversity.map <- reclassify(obj@esdms[[1]]@projection[[1]], c(-Inf,Inf, 0))

  # Useless datacheck to prevent bugs to remove after debugging
  for (i in seq_len(length(obj@esdms))) {
    if (!inherits(obj@esdms[[i]]@projection, "RasterLayer")) {
      if (verbose) {
        cat("Error", obj@esdms[[i]]@name, "is not a raster but a",
            class(obj@esdms[[i]]@projection)[1],
            ".\nIt will be removed for the stacking")
      }
      obj@esdms[[i]] <- NULL
    }
  }

  if (method == "bSSDM") {
    # Threshold and sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by thresholding and then summing. \n")
    diversity.map <- sum(stack(lapply(obj@esdms, function(x)
      reclassify(x@projection,
                 c(-Inf,x@evaluation$threshold,0,x@evaluation$threshold,Inf,1))
      )))
  }

  if (method == "pSSDM") {
    # Individual probabilities sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by summing individual probabilities. \n")
    diversity.map <- sum(stack(lapply(obj@esdms, function(x) x@projection)))
  }

  if (method == "Bernoulli") {
    # Random Bernoulli distribution (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by drawing repeatedly from a Bernoulli distribution. \n")
    proba <- stack(lapply(obj@esdms, function(x) x@projection))
    diversity.map <- calc(proba, fun = function(...) {
      x <- c(...)
      x[is.na(x)] <- 0
      return(rbinom(lengths(x), rep.B, x))
    }, forcefun = TRUE)
    diversity.map <- sum(diversity.map)/length(esdms)/rep.B
  }

  if(method=="MaximumLikelihood") {
    # Species richness adjustment by using maximum likelihood estimates for parameters on the logit-transformed occurrence probabilities (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by maximum likelihood adjustment with Calabrese bias correction. \n")
    # create empty map (keeping extent and resolution of projections)
    diversity.map <- setValues(obj@esdms[[1]]@projection,NA)
    # get observed richness (vector)
    Richness <- .richness(obj)
    ind_notNA <- Which(!is.na(Richness),cells=TRUE)
    obs.sr <- Richness[ind_notNA]
    # get predicted richness (site-by-species occurrence probabilities)
    pred.prob <- as.data.frame(sapply(obj@esdms, function(x) x@projection[ind_notNA]))
    ### Start of code from Damaris Zurell/Justin Calabrese (minor edits)
    # helper functions and maximum likelihood function
    logit = function(x) {x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x));   ;log(x/(1 - x))}
    invlogit = function(x) {exp(x)/(1+exp(x))}

    nLL.Calabrese <- function(par,sr,probs) {
      bysite <- function(j) {
        logit.probs <- logit(as.numeric(probs[j,]))
        corr.probs <- invlogit( logit.probs + par[1]*sr[j] + par[2] )
        log(dpoibin(sr[j],as.numeric(corr.probs)))
      }
      - sum(sapply(seq_len(length(sr)),bysite)) 	# optim performs minimization, so for maximum likelihood we need to invert
    }

  # find the maximum likelihood estimates for adjustment parameters
    adj.par <- optim(par=c(0,0), fn=nLL.Calabrese, sr=obs.sr, probs=pred.prob)
  # adjust the predicted probabilities, transform back with invlogit and stack
    corr.sr <- rowSums(apply(pred.prob,2,FUN=function(x){invlogit(logit(x)+adj.par$par[1]*obs.sr+adj.par$par[2])}))
    ### End of code from Damaris Zurell/Justin Calabrese

    diversity.map[ind_notNA] <- corr.sr
  }


  if (method == "PRR.MEM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from MEM. \n")
    diversity.map <- .MEM(obj,Env)@projection
    esdms <- .PRR(obj, diversity.map)
  }

  if (method == "PRR.pSSDM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from pSSDM. \n")
    diversity.map <- mapDiversity(obj, method = 'pSSDM',
                                  verbose = FALSE)$diversity.map
    esdms <- .PRR(obj, diversity.map)
  }

  return(list(
    diversity.map = diversity.map,
    esdms = esdms
  ))
})

##### Internals ####

.richness <- function(obj){
  Richness <- reclassify(obj@esdms[[1]]@projection, c(-Inf, Inf, 0))
  for (i in seq_len(length(obj@esdms)))
    Richness <- Richness + rasterize(
      sf::st_as_sf(obj@esdms[[i]]@data[1:2], coords = c("X", "Y")),
      Richness, field = obj@esdms[[i]]@data$Presence,
      background = 0)
  if (all(values(Richness) %in% c(0, 1, NA)))
    stop("Observed Richness is always equal to 1, modelled richness can't be adjusted !")
  return(Richness)
}

.MEM <- function(obj, Env){
  occ <- data.frame(rasterToPoints(.richness(obj), function(x) x > 0))
  maxOcc <- max(occ$layer) # Reucing occ for algorithms
  occ$layer <- occ$layer/max(maxOcc)
  algo <- unlist(
    strsplit(obj@esdms[[1]]@parameters$algorithms,
             ".", fixed = TRUE))[-1]
  if("MAXENT" %in% algo)
    algo <- algo[-which(algo == "MAXENT")]
  MEM <- ensemble_modelling(algorithms = algo,
    Occurrences = occ, Env = Env, Xcol = "x",
    Ycol = "y", Pcol = "layer", rep = obj@esdms[[1]]@parameters$rep,
    name = "MEM", cv = obj@esdms[[1]]@parameters$cv,
    cv.param = as.numeric(unlist(
      strsplit(obj@esdms[[1]]@parameters$cv.param,
               "|", fixed = TRUE))[-1]),
    metric = obj@esdms[[1]]@parameters$metric,
    axes.metric = obj@esdms[[1]]@parameters$axes.metric,
    ensemble.metric = unlist(
      strsplit(obj@esdms[[1]]@parameters$ensemble.metric,
               ".", fixed = TRUE))[-1],
    ensemble.thresh = as.numeric(unlist(
      strsplit(obj@esdms[[1]]@parameters$ensemble.thresh,
               "|", fixed = TRUE))[-1]),
    uncertainty = FALSE,
    weight = as.logical(obj@esdms[[1]]@parameters$weight),
    verbose = FALSE)
  MEM@projection <- MEM@projection*maxOcc
  return(MEM)
}

.PRR <- function(obj, Richness){
  # Readjust each esdm binary map
  richnesses <- values(Richness)
  names(richnesses) <- seq_len(length(richnesses))
  richnesses <- as.list(richnesses)
  probabilities <- lapply(lapply(obj@esdms, FUN = slot, name = "projection"),
                          values)
  probabilities <- lapply(probabilities, function(x) {
    names(x) <- rep(seq_len(length(probabilities[[1]])))
    return(x)
  })
  probabilities <- lapply(probabilities, `[`, names(probabilities[[1]]))
  probabilities <- apply(do.call(rbind, probabilities), 2, as.list)
  binaries <- lapply(lapply(obj@esdms, FUN = slot, name = "binary"),
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
  mapply(function(esdm, binary) {
    values(esdm@binary) <- binary
    return(esdm)
  }, esdm = obj@esdms, binary = binaries, SIMPLIFY = FALSE)
}
