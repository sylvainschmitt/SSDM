#' @include Algorithm.SDM.R ensemble.R Ensemble.SDM.R ensemble_modelling.R checkargs.R
#'   stacking.R Stacked.SDM.R
#' @importFrom shiny incProgress
#' @importFrom raster stack writeRaster
NULL

#'Build an SSDM that assembles multiple algorithms and species.
#'
#'This is a function to build an SSDM that assembles multiple algorithm and
#'species. The function takes as inputs an occurrence data frame made of
#'presence/absence or presence-only records and a raster object for data
#'extraction and projection. The function returns an S4
#'\linkS4class{Stacked.SDM} class object containing the local species richness
#'map, the between-algorithm variance map, and all evaluation tables coming with
#'(model evaluation, algorithm evaluation, algorithm correlation matrix and
#'variable importance), and a list of ensemble SDMs for each species (see
#'\code{\link{ensemble_modelling}}).
#'
#'@param algorithms character. Choice of the algorithm(s) to be run (see details
#'  below).
#'@param Occurrences data frame. Occurrence table (can be processed first by
#'  \code{\link{load_occ}}).
#'@param Env raster object. Raster object of environmental variables (can be
#'  processed first by \code{\link{load_var}}).
#'@param Xcol character. Name of the column in the occurrence table containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence. A value of 1 is presence and
#'  value of 0 is absence. If NULL presence-only dataset is assumed.
#'@param Spcol character. Name of the column containing species names or IDs.
#'@param rep integer. Number of repetitions for each algorithm.
#'@param name character. Optional name given to the final Ensemble.SDM produced.
#'@param save logical. If set to true, the SSDM is automatically saved.
#'@param path character. If save is true, the path to the directory in which the
#'  ensemble SDM will be saved.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only dataset. If PA is NULL, recommended PA selection
#'  strategy is used depending on the algorithm (see details below).
#'@param cv character. Method of cross-validation used to evaluate the ensemble
#'  SDM (see details below).
#'@param cv.param numeric. Parameters associated with the method of
#'  cross-validation used to evaluate the ensemble SDM (see details below).
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1.
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate variable relative importance (see
#'  details below).
#'@param uncertainty logical. If set to true, generates an uncertainty map and
#'  an algorithm correlation matrix.
#'@param tmp logical. If set to true, the habitat suitability map of each
#'  algorithms is saved in a temporary file to release memory. But beware: if
#'  you close R, temporary files will be deleted To avoid any loss you can
#'  save your SSDM with \code{\link{save.model}}. Depending on number,
#'  resolution and extent of models, temporary files can take a lot of disk
#'  space. Temporary files are written in R environment temporary folder.
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'@param weight logical. Choose whether or not you want the SDMs to be weighted
#'  using the selection metric or, alternatively, the mean of the selection
#'  metrics.
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random bernoulli (\strong{Bernoulli}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'@param range integer. Set a value of range restriction (in pixels) around
#'  presences occurrences on habitat suitability maps (all further points will
#'  have a null probability, see Crisp et al (2011) in references). If NULL, no
#'  range restriction will be applied.
#'@param endemism character. Define the method used to create an endemism map
#'  (see details below).
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'@param cores integer. Specify the number of CPU cores used to do the
#'  computing. You can use \code{\link[parallel]{detectCores}}) to automatically
#'  use all the available CPU cores.
#'@param ... additional parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'
#'@details \describe{ \item{algorithms}{'all' allows you to call directly all
#'  available algorithms. Currently, available algorithms include Generalized
#'  linear model (\strong{GLM}), Generalized additive model (\strong{GAM}),
#'  Multivariate adaptive regression splines (\strong{MARS}), Generalized
#'  boosted regressions model (\strong{GBM}), Classification tree analysis
#'  (\strong{CTA}), Random forest (\strong{RF}), Maximum entropy
#'  (\strong{MAXENT}), Artificial neural network (\strong{ANN}), and Support
#'  vector machines (\strong{SVM}). Each algorithm has its own parameters
#'  settable with the \strong{...} (see each algorithm section below to set
#'  their parameters).} \item{"PA"}{list with two values: \strong{nb} number of
#'  pseudo-absences selected, and \strong{strat} strategy used to select
#'  pseudo-absences: either random selection or disk selection. We set default
#'  recommendation from Barbet-Massin et al. (2012) (see reference).}
#'  \item{cv}{\strong{Cross-validation} method used to split the occurrence
#'  dataset used for evaluation: \strong{holdout} data are partitioned into a
#'  training set and an evaluation set using a fraction (\emph{cv.param[1]}) and
#'  the operation can be repeated (\emph{cv.param[2]}) times, \strong{k-fold}
#'  data are partitioned into k (\emph{cv.param[1]}) folds being k-1 times in
#'  the training set and once the evaluation set and the operation can be
#'  repeated (\emph{cv.param[2]}) times, \strong{LOO} (Leave One Out) each point
#'  is successively taken as evaluation data.} \item{metric}{Choice of the
#'  metric used to compute the binary map threshold and the confusion matrix (by
#'  default SES as recommended by Liu et al. (2005), see reference below):
#'  \strong{Kappa} maximizes the Kappa, \strong{CCR} maximizes the proportion of
#'  correctly predicted observations, \strong{TSS} (True Skill Statistic)
#'  maximizes the sum of sensitivity and specificity, \strong{SES} uses the
#'  sensitivity-specificity equality, \strong{LW} uses the lowest occurrence
#'  prediction probability, \strong{ROC} minimizes the distance between the ROC
#'  plot (receiving operating curve) and the upper left corner
#'  (1,1).}\item{axes.metric}{Choice of the metric used to evaluate the variable
#'  relative importance (difference between a full model and one with each
#'  variable successively omitted): \strong{Pearson} (computes a simple
#'  Pearson's correlation \emph{r} between predictions of the full model and the
#'  one without a variable, and returns the score \emph{1-r}: the highest the
#'  value, the more influence the variable has on the model), \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{ensemble.metric}{Ensemble metric(s) used to select SDMs: \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{method}{Choice of the method used to compute the local species
#'  richness map (see Calabrese et al. (2014) and D'Amen et al (2015) for more
#'  informations, see reference below): \strong{pSSDM} sum probabilities of
#'  habitat suitability maps, \strong{Bernoulli} drawing repeatedly from a
#'  Bernoulli distribution, \strong{bSSDM} sum the binary map obtained with the
#'  thresholding (depending on the metric, see metric parameter),
#'  \strong{MaximumLikelihood} adjust species richness using maximum likelihood parameter estimates on the logit-transformed occurrence probabilities (see Calabrese et al. (2014)), \strong{PRR.MEM} model richness with a macroecological model
#'  (MEM) and adjust each ESDM binary map by ranking habitat suitability and
#'  keeping as much as predicted richness of the MEM, \strong{PRR.pSSDM} model
#'  richness with a pSSDM and adjust each ESDM binary map by ranking habitat
#'  suitability and keeping as much as predicted richness of the pSSDM}
#'  \item{endemism}{Choice of the method used to compute the endemism map (see
#'  Crisp et al. (2001) for more information, see reference below):
#'  \strong{NULL} No endemism map, \strong{WEI} (Weighted Endemism Index)
#'  Endemism map built by counting all species in each cell and weighting each
#'  by the inverse of its range, \strong{CWEI} (Corrected Weighted Endemism
#'  Index) Endemism map built by dividing the weighted endemism index by the
#'  total count of species in the cell. First string of the character is the
#'  method either WEI or CWEI, and in those cases second string of the vector is
#'  used to precise range calculation, whether the total number of occurrences
#'  \strong{'NbOcc'} whether the surface of the binary map species distribution
#'  \strong{'Binary'}.} \item{...}{See algorithm in detail section} }
#'
#'@section Generalized linear model (\strong{GLM}) : Uses the \code{glm}
#'  function from the package 'stats', you can set the following parameters (see
#'  \code{\link[stats]{glm}} for more details): \describe{
#'  \item{test}{character. Test used to evaluate the SDM, default 'AIC'.}
#'  \item{epsilon}{numeric. Positive convergence tolerance eps ; the iterations
#'  converge when \emph{|dev - dev_{old}|/(|dev| + 0.1) < eps}. By default, set
#'  to 10e-08.} \item{maxit}{numeric. Integer giving the maximal number of IWLS
#'  (Iterative Weighted Last Squares) iterations, default 500.} }
#'
#'@section Generalized additive model (\strong{GAM}) : Uses the \code{gam}
#'  function from the package 'mgcv', you can set the following parameters (see
#'  \code{\link[mgcv]{gam}} for more details): \describe{ \item{test}{character.
#'  Test used to evaluate the model, default 'AIC'.} \item{epsilon}{numeric.
#'  This is used for judging conversion of the GLM IRLS (Iteratively Reweighted
#'  Least Squares) loop, default 10e-08.} \item{maxit}{numeric. Maximum number
#'  of IRLS iterations to perform, default 500.} }
#'
#'@section Multivariate adaptive regression splines (\strong{MARS}) : Uses the
#'  \code{earth} function from the package 'earth', you can set the following
#'  parameters (see \code{\link[earth]{earth}} for more details): \describe{
#'  \item{degree}{integer. Maximum degree of interaction (Friedman's mi) ; 1
#'  meaning build an additive model (i.e., no interaction terms). By default,
#'  set to 2.} }
#'
#'@section Generalized boosted regressions model (\strong{GBM}) : Uses the
#'  \code{gbm} function from the package 'gbm,' you can set the following
#'  parameters (see \code{\link[gbm]{gbm}} for more details): \describe{
#'  \item{trees}{integer. The total number of trees to fit. This is equivalent
#'  to the number of iterations and the number of basis functions in the
#'  additive expansion. By default, set to 2500.} \item{final.leave}{integer.
#'  minimum number of observations in the trees terminal nodes. Note that this
#'  is the actual number of observations not the total weight. By default, set
#'  to 1.} \item{algocv}{integer. Number of cross-validations, default 3.}
#'  \item{thresh.shrink}{integer. Number of cross-validation folds to perform.
#'  If cv.folds>1 then gbm, in addition to the usual fit, will perform a
#'  cross-validation. By default, set to 1e-03.} }
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the \code{rpart}
#'  function from the package 'rpart', you can set the following parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{
#'  \item{final.leave}{integer. The minimum number of observations in any
#'  terminal node, default 1.} \item{algocv}{integer. Number of
#'  cross-validations, default 3.} }
#'
#'@section Random Forest (\strong{RF}) : Uses the \code{randomForest} function
#'  from the package 'randomForest', you can set the following parameters (see
#'  \code{\link[randomForest]{randomForest}} for more details): \describe{
#'  \item{trees}{integer. Number of trees to grow. This should not be set to a
#'  too small number, to ensure that every input row gets predicted at least a
#'  few times. By default, set to 2500.} \item{final.leave}{integer. Minimum
#'  size of terminal nodes. Setting this number larger causes smaller trees to
#'  be grown (and thus take less time). By default, set to 1.} }
#'
#'@section Maximum Entropy (\strong{MAXENT}) : Uses the \code{maxent} function
#'  from the package 'dismo'. Make sure that you have correctly installed the
#'  maxent.jar file in the folder ~\\R\\library\\version\\dismo\\java available
#'  at \url{https://www.cs.princeton.edu/~schapire/maxent/} (see
#'  \code{\link[dismo]{maxent}} for more details).
#'
#'@section Artificial Neural Network (\strong{ANN}) : Uses the \code{nnet}
#'  function from the package 'nnet', you can set the following parameters (see
#'  \code{\link[nnet]{nnet}} for more details): \describe{ \item{maxit}{integer.
#'  Maximum number of iterations, default 500.} }
#'
#'@section Support vector machines (\strong{SVM}) : Uses the \code{svm} function
#'  from the package 'e1071', you can set the following parameters (see
#'  \code{\link[e1071]{svm}} for more details): \describe{ \item{epsilon}{float.
#'  Epsilon parameter in the insensitive loss function, default 1e-08.}
#'  \item{algocv}{integer. If an integer value k>0 is specified, a k-fold
#'  cross-validation on the training data is performed to assess the quality of
#'  the model: the accuracy rate for classification and the Mean Squared Error
#'  for regression. By default, set to 3.} }
#'
#'@section Warning : Depending on the raster object resolution the process can
#'  be more or less time and memory consuming.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#'
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Results plotting
#' plot(SSDM)
#' }
#'
#'@seealso \code{\link{modelling}} to build simple SDMs.
#'
#'
#'@references M. D'Amen, A. Dubuis, R. F. Fernandes, J. Pottier, L. Pelissier, &
#'  A Guisan (2015) "Using species richness and functional traits prediction to
#'  constrain assemblage predicitions from stacked species distribution models"
#'  \emph{Journal of Biogeography} 42(7):1255-1266
#'  \url{http://doc.rero.ch/record/235561/files/pel_usr.pdf}
#'
#'  M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012) "Selecting
#'  pseudo-absences for species distribution models: how, where and how many?"
#'  \emph{Methods Ecology and Evolution} 3:327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'
#'
#'  J.M. Calabrese, G. Certain, C.  Kraan, & C.F. Dormann (2014) "Stacking
#'  species distribution  models  and  adjusting  bias  by linking them to
#'  macroecological models." \emph{Global Ecology and Biogeography} 23:99-112
#'  \url{http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/calabrese2013globalecolbiogeogr.pdf}
#'
#'
#'
#'  M. D. Crisp, S. Laffan, H. P. Linder & A. Monro (2001) "Endemism in the
#'  Australian flora"  \emph{Journal of Biogeography} 28:183-198
#'  \url{http://biology-assets.anu.edu.au/hosted_sites/Crisp/pdfs/Crisp2001_endemism.pdf}
#'
#'
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting
#'  thresholds of occurrence in the prediction of species distributions."
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'
#'@export
stack_modelling <- function(algorithms,
                           # Modelling data input
                           Occurrences, Env,
                           # Occurrences reading
                           Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL, Spcol = 'SpeciesID',
                           # Model creation
                           rep = 10, name = NULL, save = FALSE, path = getwd(),
                           # Pseudo-absences definition
                           PA = NULL,
                           # Evaluation parameters
                           cv = 'holdout', cv.param = c(0.7,1), thresh = 1001,
                           axes.metric = 'Pearson', uncertainty = TRUE, tmp = FALSE,
                           # Assembling parameters
                           ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = TRUE,
                           # Diversity map computing
                           method = 'pSSDM', metric = 'SES', rep.B = 1000,
                           # Range restriction and endemism
                           range = NULL, endemism = c('WEI','Binary'),
                           # Informations parameters
                           verbose = TRUE, GUI = FALSE, cores = 1,
                           # Modelling parameters
                           ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, Spcol = Spcol, rep = rep,
             name = name, save = save, path = path, PA = PA, cv = cv, cv.param = cv.param,
             thresh = thresh, axes.metric = axes.metric, uncertainty = uncertainty,
             tmp = tmp, ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, method = method, metric = metric, rep.B = rep.B, range = range,
             endemism = endemism, verbose = verbose, GUI = GUI, cores = cores)

  # Test if algorithm is available
  available.algo <- c("GLM", "GAM", "MARS", "GBM", "CTA", "RF", "MAXENT",
                      "ANN", "SVM")
  if ("all" %in% algorithms) {
    algorithms <- available.algo
  }
  for (i in seq_len(length(algorithms))) {
    if (!(algorithms[[i]] %in% available.algo)) {
      stop(algorithms[[i]], " is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, SVM")
    }
  }
  if (tmp) {
    tmppath <- get("tmpdir", envir = .PkgEnv)
    if (!("/.esdms" %in% list.dirs(tmppath)))
      (dir.create(paste0(tmppath, "/.esdms")))
  }

  # Ensemble models creation
  if (verbose) {
    cat("#### Ensemble models creation ##### \n\n")
  }
  species <- levels(as.factor(Occurrences[, which(names(Occurrences) == Spcol)]))

  if (cores > 0 && requireNamespace("parallel", quietly = TRUE)) {
    if (verbose) {
      cat("Opening clusters,", cores, "cores \n")
    }
    if ((parallel::detectCores() - 1) < cores) {
      warning("It seems you attributed more cores than your CPU have !")
    }
    cl <- parallel::makeCluster(cores, outfile = "")
    if (verbose) {
      cat("Exporting environment to clusters \n")
    }
    parallel::clusterExport(cl, varlist = c(lsf.str(envir = globalenv()),
                                            ls(envir = environment())), envir = environment())
    esdms <- parallel::parLapply(cl, species, function(species) {
      esdm.name <- species
      Spoccurrences <- subset(Occurrences, Occurrences[which(names(Occurrences) ==
                                                               Spcol)] == species)
      if (verbose) {
        cat("Ensemble modelling :", esdm.name, "\n\n")
      }
      esdm <- try(ensemble_modelling(algorithms, Spoccurrences, Env, Xcol,
                                    Ycol, Pcol, rep = rep, name = esdm.name, save = FALSE, path = path,
                                    PA = PA, cv = cv, cv.param = cv.param, thresh = thresh, metric = metric,
                                    axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
                                    ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
                                    weight = weight, verbose = verbose, GUI = FALSE, n.cores = 1,
                                    ...))
      if (GUI) {
        incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                     Spcol)]))) + 1), detail = paste(species, " ensemble SDM built"))
      }
      if (inherits(esdm, "try-error")) {
        if (verbose) {
          cat(esdm)
        }
        esdm <- NULL
      } else {
        if (tmp && !is.null(esdm)) {
          esdm@projection <- writeRaster(esdm@projection[[1]], paste0(tmppath,
                                                                    "/.esdms/proba", esdm.name), overwrite = TRUE)
          esdm@binary <- writeRaster(esdm@binary[[1]], paste0(tmppath,
                                                            "/.esdms/bin", esdm.name), overwrite = TRUE)
          esdm@uncertainty <- writeRaster(esdm@uncertainty, paste0(tmppath,
                                                                 "/.esdms/uncert", esdm.name), overwrite = TRUE)
        }
        if (verbose) {
          cat("\n\n")
        }
      }
      return(esdm)
    })
    if (verbose) {
      cat("Closing clusters \n")
    }
    parallel::stopCluster(cl)

  } else {
    esdms <- lapply(species, function(species) {
      esdm.name <- species
      Spoccurrences <- subset(Occurrences, Occurrences[which(names(Occurrences) ==
                                                               Spcol)] == species)
      if (verbose) {
        cat("Ensemble modelling :", esdm.name, "\n\n")
      }
      esdm <- try(ensemble_modelling(algorithms, Spoccurrences, Env, Xcol,
                                    Ycol, Pcol, rep = rep, name = esdm.name, save = FALSE, path = path,
                                    PA = PA, cv = cv, cv.param = cv.param, thresh = thresh, metric = metric,
                                    axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
                                    ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
                                    weight = weight, verbose = verbose, GUI = FALSE, ...))
      if (GUI) {
        incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                     Spcol)]))) + 1), detail = paste(species, " ensemble SDM built"))
      }
      if (inherits(esdm, "try-error")) {
        if (verbose) {
          cat(esdm)
        }
        esdm <- NULL
      } else {
        if (tmp && !is.null(esdm)) {
          esdm@projection <- writeRaster(esdm@projection[[1]], paste0(tmppath,
                                                                    "/.esdms/proba", esdm.name), overwrite = TRUE)
          esdm@uncertainty <- writeRaster(esdm@uncertainty, paste0(tmppath,
                                                                 "/.esdms/uncert", esdm.name), overwrite = TRUE)
        }
        if (verbose) {
          cat("\n\n")
        }
      }
      return(esdm)
    })
  }

  esdms <- esdms[!sapply(esdms, is.null)]

  # Species stacking
  if (length(esdms) < 2) {
    if (verbose) {
      stop("Less than two species models were retained, you should lower the ensemble threshold value (ensemble.thresh parameter).")
    } else {
      return(NULL)
    }
  } else {
    if (verbose) {
      cat("#### Species stacking with ensemble models ##### \n\n")
    }
    if (!is.null(name)) {
      esdms["name"] <- name
    }
    esdms["method"] <- method
    esdms["rep.B"] <- rep.B
    if (method %in% c("PRR.MEM", "PRR.pSSDM")) {
      esdms["Env"] <- Env
    }
    if (!is.null(range)) {
      esdms["range"] <- range
    }
    esdms$endemism <- endemism
    esdms["verbose"] <- verbose
    stack <- do.call(stacking, esdms)
  }

  if (!is.null(stack)) {
    # Paremeters
    stack@parameters$sp.nb.origin <- length(levels(as.factor(Occurrences[,
                                                                         which(names(Occurrences) == Spcol)])))
    if (GUI) {
      incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                   Spcol)]))) + 1), detail = "SSDM built")
    }

    # Saving
    if (save) {
      if (verbose) {
        cat("#### Saving ##### \n\n")
      }
      if (!is.null(name)) {
        save.stack(stack, name = name, path = path)
      } else {
        save.stack(stack, path = path)
      }
    }
  }
  return(stack)
}
