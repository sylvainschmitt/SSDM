#' @include Ensemble.SDM.R checkargs.R
#' @importFrom sf st_as_sf st_buffer
#' @importFrom raster raster stack reclassify mask calc overlay values rasterize rasterToPoints values<-
#' @importFrom stats lm
NULL

#'Stack different ensemble SDMs in an SSDM
#'
#'This is a function to stack several ensemble SDMs in an SSDM. The function
#'takes as inputs several S4 \linkS4class{Ensemble.SDM} class objects produced
#'with \code{\link{ensemble_modelling}} or \code{\link{ensemble}} functions. The
#'function returns an S4 \linkS4class{Stacked.SDM} class object containing the
#'local species richness map, the between-algorithm variance map, and all
#'evaluation tables coming with (model evaluation, algorithm evaluation,
#'algorithm correlation matrix and variable importance), and a list of ensemble
#'SDMs for each species (see \code{\link{ensemble_modelling}}).
#'
#'@param esdm,... character. Ensemble SDMs to be stacked.
#'@param name character. Optional name given to the final SSDM produced (by
#'  default 'Species.SDM').
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random bernoulli (\strong{Bernoulli}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'@param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}). Needed only for stacking
#'  method using probability ranking from richness (\strong{PRR}).
#'@param range integer. Set a value of range restriction (in pixels) around
#'  presences occurrences on habitat suitability maps (all further points will
#'  have a null probability, see Crisp et al (2011) in references). If NULL, no
#'  range restriction will be applied.
#'@param endemism character. Define the method used to create an endemism map
#'  (see details below).
#'@param eval logical. If set to FALSE, disable stack evaluation.
#' @param uncertainty logical. If set to TRUE, generates an uncertainty map and
#'  an algorithm correlation matrix.
#'@param verbose logical. If set to TRUE, allows the function to print text in
#'  the console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
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
#'  \strong{Endemism:} Choice of the method used to compute the endemism map
#'  (see Crisp et al. (2001) for more information, see reference below):
#'  \describe{\item{NULL}{No endemism map}\item{WEI}{(Weighted Endemism Index)
#'  Endemism map built by counting all species in each cell and weighting each
#'  by the inverse of its range} \item{CWEI}{(Corrected Weighted Endemism Index)
#'  Endemism map built by dividing the weighted endemism index by the total
#'  count of species in the cell.}}First string of the character is the method
#'  either WEI or CWEI, and in those cases second string of the vector is used
#'  to precise range calculation, whether the total number of occurrences
#'  \strong{'NbOcc'} whether the surface of the binary map species distribution
#'  \strong{'Binary'}.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occ1 <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#' Occ2 <- subset(Occurrences, Occurrences$SPECIES == 'gracilis')
#'
#' # SSDM building
#' ESDM1 <- ensemble_modelling(c('CTA', 'SVM'), Occ1, Env, rep = 1,
#'                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                            name = 'elliptica', ensemble.thresh = c(0.6))
#' ESDM2 <- ensemble_modelling(c('CTA', 'SVM'), Occ2, Env, rep = 1,
#'                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                            name = 'gracilis', ensemble.thresh = c(0.6))
#' SSDM <- stacking(ESDM1, ESDM2)
#'
#' # Results plotting
#' plot(SSDM)
#' }
#'
#'@seealso \code{\link{stack_modelling}} to build SSDMs.
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
#'  M. D. Crisp, S. Laffan, H. P. Linder & A. Monro (2001) "Endemism in the
#'  Australian flora"  \emph{Journal of Biogeography} 28:183-198
#'  \url{http://biology-assets.anu.edu.au/hosted_sites/Crisp/pdfs/Crisp2001_endemism.pdf}
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting
#'  thresholds of occurrence in the prediction of species distributions."
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'@rdname stacking
#'@export
setGeneric('stacking', function(esdm, ..., name = NULL, method = 'pSSDM', rep.B = 1000,
                                Env = NULL, range = NULL, endemism = c('WEI','Binary'), eval = TRUE, uncertainty=TRUE,
                                verbose = TRUE, GUI = FALSE) {return(standardGeneric('stacking'))})

#' @rdname stacking
#' @export
setMethod('stacking', 'Ensemble.SDM', function(esdm, ..., name = NULL, method = 'pSSDM', rep.B = 1000,
                                               Env = NULL, range = NULL, endemism = c('WEI','Binary'), eval = TRUE, uncertainty=TRUE,
                                               verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(esdm = esdm, name = name, method = method, rep.B = rep.B, range = range,
             endemism = endemism, eval = eval, verbose = verbose, GUI = GUI)

  esdms <- list(esdm, ...)
  if (length(esdms) < 2) {
    stop("You neeed more than one ensemble SDM to do stackings")
  }
  names <- c()
  for (i in seq_len(length(esdms))) {
    if (esdms[[i]]@name %in% names) {
      stop("Ensemble models can't have the same name, you need to rename one of ",
           esdms[[i]]@name)
    } else {
      names <- c(names, esdms[[i]]@name)
    }
  }
  if (verbose) {
    cat("Stack creation... \n")
  }
  stack <- Stacked.SDM(diversity.map = reclassify(esdm@projection[[1]], c(-Inf,Inf, 0)),
                       endemism.map = reclassify(esdm@projection[[1]], c(-Inf, Inf, 0)),
                       uncertainty = if(!uncertainty){raster()}else{reclassify(esdm@uncertainty, c(-Inf, Inf, NA))},
                       parameters = esdm@parameters)

  # ESDMs
  for (i in seq_len(length(esdms))) {
    suppressWarnings({
      stack@esdms[esdms[[i]]@name] <- esdms[[i]]
    })
  }

  # Name
  if (verbose) {
    cat("   naming...")
  }
  if (is.null(name)) {
    name <- "Species"
  }
  stack@name <- paste0(name, ".SSDM")
  if (verbose) {
    cat(" done. \n")
  }

  # Range restriction
  if (verbose) {
    cat("   range restriction...")
  }
  if (!is.null(range)) {
    for (j in seq_len(length(esdms))) {
      nbocc <- length(as.factor(esdms[[j]]@data$Presence[esdms[[j]]@data$Presence ==
                                                          1]))/sum(esdms[[j]]@algorithm.evaluation$kept.model)
      occ <- esdms[[j]]@data[1:nbocc, ]
      occ <- occ[which(occ$Presence == 1), 1:2]
      circles <- st_as_sf(occ, coords = c("X", "Y"))
      circles <- st_buffer(circles, range * res)
      esdms[[j]]@projection <- mask(esdms[[j]]@projection, circles, updatevalue = 0)
    }
  }
  if (verbose) {
    cat(" done. \n")
  }

  # Diversity map
  if (verbose)
    cat("   diversity mapping...")
  diversity <- mapDiversity(stack, method, rep.B, verbose, Env)
  stack@diversity.map <- diversity$diversity.map
  if(!is.null(diversity$esdms))
    stack@esdms <- diversity$esdms
  names(stack@diversity.map) <- "diversity"
  if (verbose)
    cat(" done. \n")

  # uncertainty map
  if(uncertainty){
    if (verbose) {
      cat("   uncertainty mapping...")
    }
    uncertainties <- stack()
    for (i in seq_len(length(esdms))) {
      a <- try(esdms[[i]]@uncertainty)
      if (inherits(a, "try-error")) {
        if (verbose) {
          cat("Ensemble model", esdms[[i]]@name, "uncertainty map not computed")
        }
      } else {
        b <- try(stack(uncertainties, a))
        if (inherits(b, "try-error")) {
          if (verbose) {
            cat("Ensemble model", esdms[[i]]@name, ":", b)
          }
        } else {
          uncertainties <- b
        }
      }
    }
    a <- try(calc(uncertainties, mean))
    if (inherits(a, "try-error")) {
      if (verbose) {
        cat("No uncertainty map to do uncertainty mapping")
      }
    } else {
      stack@uncertainty <- a
      names(stack@uncertainty) <- "uncertainty"
    }
    if (verbose) {
      cat(" done. \n")
    }
  }

  # endemism map
  if (verbose) {
    cat("   endemism mapping...")
  }
  if (is.null(endemism)) {
    "deactivated"
  } else {
    for (i in seq_len(length(esdms))) {
      if (endemism[2] == "NbOcc") {
        endweight <- length(as.factor(esdms[[i]]@data$Presence[esdms[[i]]@data$Presence ==
                                                                1]))/sum(esdms[[i]]@algorithm.evaluation$kept.model)
      } else if (endemism[2] == "Binary") {
        endweight <- sum(values(reclassify(esdms[[i]]@projection, c(-Inf,
                                                                   esdms[[i]]@evaluation$threshold, 0, esdms[[i]]@evaluation$threshold,
                                                                   Inf, 1))), na.rm = TRUE)
      }
      if (endemism[1] == "WEI") {
        stack@endemism.map <- stack@endemism.map + esdms[[i]]@projection/endweight
      } else if (endemism[1] == "CWEI") {
        stack@endemism.map <- stack@endemism.map + overlay(esdms[[i]]@projection,
                                                           stack@diversity.map, fun = function(x, y) {
                                                             y <- round(y)
                                                             x[which(y > 0)] <- x[which(y > 0)]/endweight/y[which(y >
                                                                                                                    0)]
                                                             return(x)
                                                           })
      }
    }
    stack@endemism.map <- stack@endemism.map/stack@endemism.map@data@max
  }
  if (verbose) {
    cat(" done. \n")
  }

  # variable Importance
  if (verbose) {
    cat("   comparing variable importance...")
  }
  stack@variable.importance <- esdm@variable.importance
  for (i in 2:length(esdms)) {
    a <- try(rbind(stack@variable.importance, esdms[[i]]@variable.importance))
    if (inherits(a, "try-error")) {
      cat(a)
    } else {
      stack@variable.importance <- a
    }
  }
  a <- stack@variable.importance[1:2, ]
  row.names(a) <- c("Mean", "SD")
  for (i in seq_len(length(stack@variable.importance))) {
    a[i] <- c(mean(stack@variable.importance[, i]), sd(stack@variable.importance[,
                                                                                 i]))
  }
  stack@variable.importance <- a
  if (verbose) {
    cat(" done. \n")
  }

  # Algorithm Correlation
  if (verbose) {
    cat("   comparing algorithms correlation...")
  }
  algo <- c()  # Listing all algorithms presents in esdms and renaming esdms row and columns
  for (i in seq_len(length(esdms))) {
    if (length(esdms[[i]]@algorithm.correlation) == 0) {
      if (verbose) {
        cat("\n", esdms[[i]]@name, "algorithms correlation has not been computed. \n")
      }
    } else {
      for (j in seq_len(length(esdms[[i]]@algorithm.correlation))) {
        if (length(strsplit(names(esdms[[i]]@algorithm.correlation)[j],
                            ".", fixed = TRUE)[[1]]) > 1) {
          a <- strsplit(row.names(esdms[[i]]@algorithm.correlation)[j],
                        ".SDM", fixed = TRUE)[[1]][1]
          a <- tail(strsplit(a, ".", fixed = TRUE)[[1]], n = 1)
          names(esdms[[i]]@algorithm.correlation)[j] <- a
          row.names(esdms[[i]]@algorithm.correlation)[j] <- a
        }
        if (!(names(esdms[[i]]@algorithm.correlation)[j] %in% algo)) {
          algo <- c(algo, names(esdms[[i]]@algorithm.correlation)[j])
        }
      }
    }
  }
  mcorr <- data.frame(matrix(nrow = length(algo), ncol = length(algo)))
  names(mcorr) <- algo
  row.names(mcorr) <- algo
  if (length(algo) > 0) {
    for (i in seq_len(length(algo))) {
      for (j in seq_len(length(algo))) {
        if (i > j) {
          corr <- c()
          for (k in seq_len(length(esdms))) {
            if (length(esdms[[k]]@algorithm.correlation) != 0) {
              row <- which(row.names(esdms[[k]]@algorithm.correlation) ==
                             row.names(mcorr)[j])
              col <- which(names(esdms[[k]]@algorithm.correlation) ==
                             names(mcorr)[i])
              if (length(row) > 0 && length(col) > 0) {
                corr <- c(corr, esdms[[k]]@algorithm.correlation[row,
                                                                col])
              }
            }
            mcorr[i, j] <- mean(corr, na.rm = TRUE)
          }
        }
      }
    }
  }
  stack@algorithm.correlation <- mcorr
  if (verbose) {
    cat(" done. \n")
  }

  # Algorithm Evaluation
  if (verbose) {
    cat("   comparing algorithms evaluation")
  }
  stack@algorithm.evaluation <- esdm@algorithm.evaluation
  for (i in 2:length(esdms)) {
    stack@algorithm.evaluation <- rbind(stack@algorithm.evaluation, esdms[[i]]@algorithm.evaluation)
  }
  # stack@algorithm.evaluation$algo <- "algo"
  stack@algorithm.evaluation$algo <- unlist(regmatches(row.names(stack@algorithm.evaluation), gregexpr("ANN|CTA|GAM|GBM|GLM|MARS|MAXENT|RF|SVM", row.names(stack@algorithm.evaluation))))
  # for (i in seq_len(length(row.names(stack@algorithm.evaluation)))) {
  #   stack@algorithm.evaluation$algo[i] <- strsplit(row.names(stack@algorithm.evaluation),
  #                                                  ".", fixed = TRUE)[[i]][2]
  # }
  stack@algorithm.evaluation <- aggregate.data.frame(stack@algorithm.evaluation[-which(names(stack@algorithm.evaluation) ==
                                                                                         "algo")],
                                                     by = list(stack@algorithm.evaluation[, which(names(stack@algorithm.evaluation) ==
                                                                                                    "algo")]), FUN = mean)
  row.names(stack@algorithm.evaluation) <- stack@algorithm.evaluation$Group.1
  stack@algorithm.evaluation <- stack@algorithm.evaluation[-1]

  if (verbose) {
    cat(" done. \n")
  }

  # Evaluation
  if (eval){
    if (verbose)
      cat("   evaluating...")
    stack@evaluation <- evaluate(stack)
  }

  # Parameters
  stack@parameters$method <- method
  if (method == "B") {
    stack@parameters$rep.B <- rep.B
  }
  stack@parameters$range <- range
  stack@parameters$endemism <- paste0(endemism[1], "|", endemism[2])

  return(stack)
})
