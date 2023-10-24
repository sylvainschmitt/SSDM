#' @include modelling.R ensemble.R Ensemble.SDM.R checkargs.R
#' @importFrom shiny incProgress
#' @importFrom raster writeRaster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators icount
#' @importFrom parallel detectCores makeCluster stopCluster
NULL

#'Build an ensemble SDM that assembles multiple algorithms
#'
#'Build an ensemble SDM that assembles multiple algorithms for a single species.
#'The function takes as inputs an occurrence data frame made of presence/absence
#'or presence-only records and a raster object for data extraction and
#'projection. The function returns an S4 \linkS4class{Ensemble.SDM} class object
#'containing the habitat suitability map, the binary map, the between-algorithm
#'variance map and the associated evaluation tables (model evaluation, algorithm
#'evaluation, algorithm correlation matrix and variable importance).
#'
#'@param algorithms character. A character vector specifying the algorithm
#'  name(s) to be run (see details below).
#'@param Occurrences data frame. Occurrences table (can be processed first by
#'  \code{\link{load_occ}}).
#'@param Env raster object. RasterStack object of environmental variables (can
#'  be processed first by \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurrence table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence. A value of 1 is presence and
#'  value of 0 is absence. If NULL presence-only dataset is assumed.
#'@param rep integer. Number of repetitions for each algorithm.
#'@param name character. Optional name given to the final Ensemble.SDM produced
#'  (by default 'Ensemble.SDM').
#'@param save logical. If \code{TRUE}, the ensemble SDM is automatically saved.
#'@param path character. If save is If \code{TRUE}, the path to the directory in
#'  which the ensemble SDM will be saved.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only dataset. If PA is NULL, recommended PA selection
#'  strategy is used depending on the algorithm (see details below).
#'@param cv character. Method of cross-validation used to evaluate the ensemble
#'  SDM (see details below).
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the ensemble SDM (see details below).
#'@param final.fit.data strategy used for fitting the final/evaluated Algorithm.SDMs: 'holdout'= use same train and test data as in (last) evaluation, 'all'= train model with all data (i.e. no test data) or numeric (0-1)= sample a custom training fraction (left out fraction is set aside as test data)
#' @param bin.thresh character. Classification threshold (\code{\link[dismo]{threshold}}) used to binarize model predictions into presence/absence and compute the confusion matrix (see details below).
#' @param metric (deprecated) character. Classification threshold (\code{SDMTools::optim.thresh}) used to binarize model predictions into presence/absence and compute the confusion matrix (see details below). This argument is only kept for backwards compatibility, if possible please use \code{bin.thresh} instead.
#' @param thresh (deprecated) integer. Number of equally spaced thresholds in the interval 0-1 (\code{SDMTools::optim.thresh}). Only needed when \code{metric} is set.
#'@param axes.metric Metric used to evaluate variable relative importance (see
#'  details below).
#'@param uncertainty logical. If \code{TRUE}, generates an uncertainty map and
#'  an algorithm correlation matrix.
#'@param SDM.projections logical. If FALSE (default), the Algorithm.SDMs inside the 'sdms' slot will not contain projections (for memory saving purposes).
#'@param tmp logical or character. If \code{FALSE}, no temporary rasters are written (this could quickly fill up your working memory, if many replicates are modelled). If \code{TRUE}, temporary rasters are written to the „tmp“ directory of your R environment. If \code{character}, temporary rasters are written to a custom path. But beware: if you
#'  close R, temporary files will be deleted. To avoid any loss you can save your
#'  ensemble SDM with \code{\link{save.model}}. Depending on number, resolution
#'  and extent of models, temporary files can take a lot of disk space.
#'  Temporary files are written to the R environment temporary folder.
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'@param weight logical. If \code{TRUE}, SDMs are weighted using the ensemble
#'  metric or, alternatively, the mean of the selection metrics.
#'@param cores integer. Specify the number of CPU cores used to do the
#'  computing. You can use \code{\link[parallel]{detectCores}}) to automatically
#'@param parmode character. Parallelization mode: along 'algorithms' or 'replicates'. Defaults to 'replicates'.
#'@param verbose logical. If \code{TRUE}, allows the function to print text in
#'  the console.
#'@param GUI logical. Do not take this argument into account (parameter for the
#'  user interface).
#'@param ... additional parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Ensemble.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'
#'@details \describe{ \item{algorithms}{'all' calls all the following
#'  algorithms. Algorithms include Generalized linear model (\strong{GLM}),
#'  Generalized additive model (\strong{GAM}), Multivariate adaptive regression
#'  splines (\strong{MARS}), Generalized boosted regressions model
#'  (\strong{GBM}), Classification tree analysis (\strong{CTA}), Random forest
#'  (\strong{RF}), Maximum entropy (\strong{MAXENT}), Artificial neural network
#'  (\strong{ANN}), and Support vector machines (\strong{SVM}). Each algorithm
#'  has its own parameters settable with the \strong{...} (see each algorithm
#'  section below to set their parameters).} \item{"PA"}{list with two values:
#'  \strong{nb} number of pseudo-absences selected, and \strong{strat} strategy
#'  used to select pseudo-absences: either random selection or disk selection.
#'  We set default recommendation from Barbet-Massin et al. (2012) (see
#'  reference).} \item{cv}{\strong{Cross-validation} method used to split the
#'  occurrence dataset used for evaluation: \strong{holdout} data are
#'  partitioned into a training set and an evaluation set using a fraction
#'  (\emph{cv.param[1]}) and the operation can be repeated (\emph{cv.param[2]})
#'  times, \strong{k-fold} data are partitioned into k (\emph{cv.param[1]})
#'  folds being k-1 times in the training set and once the evaluation set and
#'  the operation can be repeated (\emph{cv.param[2]}) times, \strong{LOO}
#'  (Leave One Out) each point is successively taken as evaluation data.}
#'  \item{metric}{Choice of the metric used to compute the binary map threshold
#'  and the confusion matrix (by default SES as recommended by Liu et al.
#'  (2005), see reference below): \strong{Kappa} maximizes the Kappa,
#'  \strong{CCR} maximizes the proportion of correctly predicted observations,
#'  \strong{TSS} (True Skill Statistic) maximizes the sum of sensitivity and
#'  specificity, \strong{SES} uses the sensitivity-specificity equality,
#'  \strong{LW} uses the lowest occurrence prediction probability, \strong{ROC}
#'  minimizes the distance between the ROC plot (receiving operating
#'  characteristic curve) and the upper left corner
#'  (1,1).}\item{axes.metric}{Metric used to evaluate the variable relative
#'  importance (difference between a full model and one with each variable
#'  successively omitted): \strong{Pearson} (computes a simple Pearson's
#'  correlation \emph{r} between predictions of the full model and the one
#'  without a variable, and returns the score \emph{1-r}: the highest the value,
#'  the more influence the variable has on the model), \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{ensemble.metric}{Ensemble metric(s) used to select SDMs: \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{"..."}{See algorithm in detail section} }
#'
#'@section Generalized linear model (\strong{GLM}) : Uses the \code{glm}
#'  function from the package 'stats'. You can set parameters by supplying \code{glm.args=list(arg1=val1,arg2=val2)} (see \code{\link[stats]{glm}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{test}{character. Test used to evaluate the SDM, default 'AIC'.}
#'  \item{control}{list (created with \code{\link[stats]{glm.control}}).
#'  Contains parameters for controlling the fitting process. Default is \code{glm.control(epsilon = 1e-08, maxit = 500)}.
#'  'epsilon' is a numeric and defines the positive convergence tolerance (eps). The iterations converge when \emph{|dev - dev_old|/(|dev| + 0.1) < eps}.
#'  'maxit' is an integer giving the maximal number of IWLS (Iterative Weighted Last Squares) iterations.} }
#'
#'@section Generalized additive model (\strong{GAM}) : Uses the \code{gam}
#'  function from the package 'mgcv'. You can set parameters by supplying \code{gam.args=list(arg1=val1,arg2=val2)} (see \code{\link[mgcv]{gam}} for all settable arguments).
#'  The following parameters have defaults: \describe{\item{test}{character.
#'  Test used to evaluate the model, default 'AIC'.} \item{control}{list (created with \code{\link[mgcv]{gam.control}}).
#'  Contains parameters for controlling the fitting process. Default is \code{gam.control(epsilon = 1e-08, maxit = 500)}.
#'  'epsilon' is a numeric used for judging the conversion of the GLM IRLS (Iteratively Reweighted Least Squares) loop. 'maxit' is an integer giving the maximum number of IRLS iterations to perform.} }
#'
#'@section Multivariate adaptive regression splines (\strong{MARS}) : Uses the
#'  \code{earth} function from the package 'earth'. You can set parameters by supplying \code{mars.args=list(arg1=val1,arg2=val2)} (see \code{\link[earth]{earth}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{degree}{integer. Maximum degree of interaction (Friedman's mi) ; 1
#'  meaning build an additive model (i.e., no interaction terms). By default,
#'  set to 2.} }
#'
#'@section Generalized boosted regressions model (\strong{GBM}) : Uses the
#'  \code{gbm} function from the package 'gbm'. You can set parameters by supplying \code{gbm.args=list(arg1=val1,arg2=val2)} (see \code{\link[gbm]{gbm}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{distribution}{character. Automatically detected from the format of the presence column in the occurrence dataset.}
#'  \item{n.trees}{integer. The total number of trees to fit. This is equivalent
#'  to the number of iterations and the number of basis functions in the
#'  additive expansion. By default, set to 2500.}
#'  \item{n.minobsinnode}{integer.
#'  minimum number of observations in the trees terminal nodes. Note that this
#'  is the actual number of observations, not the total weight. By default, set
#'  to 1.}
#'  \item{cv.folds}{integer. Number of cross-validation folds to perform.
#'  If cv.folds>1 then gbm - in addition to the usual fit - will perform a
#'  cross-validation. By default, set to 3.}
#'  \item{shrinkage}{numeric. A shrinkage parameter applied to each tree in the expansion (also known as learning rate or step-size reduction). By default, set to 0.001.}
#'  \item{bag.fraction}{numeric. Fraction of the training set observations randomly selected to propose the next tree in the expansion.}
#'  \item{train.fraction}{numeric. Training fraction used to fit the first gbm. The remainder is used to compute out-of-sample estimates of the loss function. By default, set to 1 (since evaluation/holdout is done with \code{SSDM::evaluate}.}
#'  \item{n.cores}{integer. Number of cores to use for parallel computation of the CV folds. By default, set to 1. If you intend to use this, please set \code{ncores=0} to avoid conflicts.} }
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the \code{rpart}
#'  function from the package 'rpart'. You can set parameters by supplying \code{cta.args=list(arg1=val1,arg2=val2)} (see \code{\link[rpart]{rpart}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{control}{list (created with \code{\link[rpart]{rpart.control}}).
#'  Contains parameters for controlling the rpart fit. The default is \code{rpart.control(minbucket=1, xval=3)}.
#'  'mibucket' is an integer giving the minimum number of observations in any
#'  terminal node. 'xval' is an integer defining the number of
#'  cross-validations.} }
#'
#'@section Random Forest (\strong{RF}) : Uses the \code{randomForest} function
#'  from the package 'randomForest'. You can set parameters by supplying \code{cta.args=list(arg1=val1,arg2=val2)} (see \code{\link[randomForest]{randomForest}} all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{ntree}{integer. Number of trees to grow. This should not be set to a
#'  too small number, to ensure that every input row gets predicted at least a
#'  few times. By default, set to 2500.}
#'  \item{nodesize}{integer. Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). By default, set to 1.} }
#'
#'@section Maximum Entropy (\strong{MAXENT}) : Uses the \code{maxent} function
#'  from the package 'dismo'. Make sure that you have correctly installed the
#'  maxent.jar file in the folder ~\\R\\library\\version\\dismo\\java available
#'  at \url{https://biodiversityinformatics.amnh.org/open_source/maxent/}. As with the other algorithms, you can set parameters by supplying \code{maxent.args=list(arg1=val1,arg2=val2)}. Mind that arguments are passed from dismo to the MAXENT software again as an argument list  (see \code{\link[dismo]{maxent}} for more details).
#'  No specific defaults are set with this method.
#'
#'@section Artificial Neural Network (\strong{ANN}) : Uses the \code{nnet}
#'  function from the package 'nnet'. You can set parameters by supplying \code{ann.args=list(arg1=val1,arg2=val2)} (see \code{\link[nnet]{nnet}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{size}{integer. Number of units in the hidden layer. By default, set to 6.}
#'  \item{maxit}{integer. Maximum number of iterations, default 500.} }
#'
#'@section Support vector machines (\strong{SVM}) : Uses the \code{svm} function
#'  from the package 'e1071'. You can set parameters by supplying \code{svm.args=list(arg1=val1,arg2=val2)} (see \code{\link[e1071]{svm}} for all settable arguments).
#'  The following parameters have defaults: \describe{
#'  \item{type}{character. Regression/classification type SVM should be used with. By default, set to "eps-regression".}
#'  \item{epsilon}{float. Epsilon parameter in the insensitive loss function, default 1e-08.}
#'  \item{cross}{integer. If an integer value k>0 is specified, a k-fold
#'  cross-validation on the training data is performed to assess the quality of
#'  the model: the accuracy rate for classification and the Mean Squared Error
#'  for regression. By default, set to 3.}
#'  \item{kernel}{character. The kernel used in training and predicting. By default, set to "radial".}
#'  \item{gamma}{numeric. Parameter needed for all kernels, default \code{1/(length(data) -1)}.} }
#'
#'@section Warning : Depending on the raster object resolution the process can
#'  be more or less time and memory consuming.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#'
#' # ensemble SDM building
#' ESDM <- ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 1,
#'                           Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                           ensemble.thresh = c(0.6))
#'
#' # Results plotting
#' plot(ESDM)
#' }
#'
#'@seealso \code{\link{modelling}} to build SDMs with a single algorithm,
#'  \code{\link{stack_modelling}} to build SSDMs.
#'
#'@references M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012)
#'  "Selecting pseudo-absences for species distribution models: how, where and
#'  how many?" \emph{Methods Ecology and Evolution} 3:327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting
#'  thresholds of occurrence in the prediction of species distributions."
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@export
ensemble_modelling <- function(algorithms,
                              # Modelling data input
                              Occurrences, Env,
                              # Occurrences reading
                              Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                              # Model creation
                              rep = 10, name = NULL, save = FALSE, path = getwd(), cores=0, parmode = 'replicates',
                              # Pseudo-absences definition
                              PA = NULL,
                              # Evaluation parameters
                              cv = 'holdout', cv.param = c(0.7,1), final.fit.data='all', bin.thresh='SES', metric = NULL, thresh = 1001,
                              axes.metric = 'Pearson', uncertainty = TRUE, tmp = FALSE, SDM.projections=FALSE,
                              # Assembling parameters
                              ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = TRUE,
                              # Informations parameters
                              verbose = TRUE, GUI = FALSE,
                              # Modelling parameters
                              ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, rep = rep, name = name, save = save,
             path = path, PA = PA, cv = cv, cv.param = cv.param, final.fit.data = final.fit.data, bin.thresh = bin.thresh, metric = metric, thresh = thresh, axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
             ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, verbose = verbose, GUI = GUI)

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
  if (isTRUE(tmp)) {
    tmppath <- get("tmpdir", envir = .PkgEnv)
    if (!("/.models" %in% list.dirs(path)))
      (dir.create(paste0(tmppath, "/.models")))
    tmppath <- paste0(tmppath, "/.models")
  }
  if(is.character(tmp)){
    if (!("/.models" %in% tmp))
      (dir.create(paste0(tmp, "/.models")))
    tmppath <- paste0(tmp, "/.models")
  }

  # Algorithms models creation
  if (is.null(name)) {
    spname <- "species"
  } else {
    spname <- name
  }
  if (verbose) {
    cat(sprintf("#### Algorithms models creation for %s ##### %s \n\n",
                spname, format(Sys.time(), "%Y-%m-%d %T")))
  }

  models <- list()

  if(cores > 0) {
    if ((detectCores() - 1) < cores) {
      cores <- detectCores()-1
      warning(paste("It seems you attributed more cores than your CPU has! Automatic reduction to",
                    cores, "cores."))
    }
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    if(verbose){
      cat("Opening clusters,", cores, "cores \n")
    }

    ### REPLICATES PARALLELIZATION MODE
    if(parmode=="replicates"){
      for (i in seq_len(length(algorithms))){
        model.name <- paste0(algorithms[i])
        if (verbose) {
          cat("Modelling :", model.name, "\n\n")
        }
        modelrep <- foreach::foreach(iterators::icount(rep),
                                     .packages = c("raster","SSDM"),
                                     .verbose=verbose) %dopar% {
          model <- try(modelling(algorithms[i], Occurrences, Env, Xcol = Xcol,
                                 Ycol = Ycol, Pcol = Pcol, name = NULL, PA = PA,
                                 cv = cv, cv.param = cv.param, final.fit.data = final.fit.data,
                                 bin.thresh = bin.thresh, metric = metric, thresh = thresh,
                                 axes.metric = axes.metric,
                                 select = FALSE, select.metric = ensemble.metric,
                                 select.thresh = ensemble.thresh,
                                 verbose = verbose, GUI = GUI))

          if (inherits(model, "try-error")) {
            if (verbose) {
              cat(model)
            }
          } else {
            ### needs further testing
            if (!isFALSE(tmp) & !is.null(model)) {
              model@projection <- writeRaster(model@projection, paste0(tmppath,
                                                                       "/proba",
                                                                       model.name,
                                                                       Sys.getpid(),
                                                                       gsub(" |:|-","", Sys.time())),
                                              format='raster',overwrite = TRUE)
              model@binary <- writeRaster(model@binary, paste0(tmppath,
                                                               "/bin",
                                                               model.name,
                                                               Sys.getpid(),
                                                               gsub(" |:|-","", Sys.time())),
                                          format='raster', overwrite = TRUE)
            }
          }
          return(model)
        }
        suppressWarnings({
          models <- c(models,modelrep)
        })
        if (verbose) {
          cat(sprintf("%s done for %s %s \n\n", model.name, spname, format(Sys.time(), "%Y-%m-%d %T")))
        }
        if (GUI) {
          incProgress(1/(length(algorithms) + 1), detail = paste(algorithms[i],
                                                                 "SDM built"))
        }
      }
    } # end parmode replicates

    ### ALGORITHM PARALLELIZATION MODE
    if(parmode=="algorithms"){
        if (verbose) {
          cat(paste(Sys.time(),"Start modelling", paste0(algorithms,collapse="/"), "in parallel \n\n"))
        }
        models <- foreach::foreach(i=1:length(algorithms),.packages = c("raster","SSDM"),
                                   .verbose=verbose) %dopar% {
          modelrep <- list()
          for (j in 1:rep) {
            model.name <- paste0(algorithms[i], ".", sprintf(paste0("%0",nchar(rep),"d"),j))
            model <- try(modelling(algorithms[i], Occurrences, Env, Xcol = Xcol,
                                 Ycol = Ycol, Pcol = Pcol, name = NULL, PA = PA, cv = cv,
                                 cv.param = cv.param, final.fit.data = final.fit.data,
                                 bin.thresh = bin.thresh, metric = metric, thresh = thresh,
                                 axes.metric = axes.metric,
                                 select = FALSE, select.metric = ensemble.metric,
                                 select.thresh = ensemble.thresh,
                                 verbose = verbose, GUI = GUI))

          if (inherits(model, "try-error")) {
            if (verbose) {
              cat(model)
            }
          } else {
            ### needs further testing
            if (!isFALSE(tmp) & !is.null(model)) {
              model@projection <- writeRaster(model@projection, paste0(tmppath,
                                                                       "/proba", model.name,
                                                                       Sys.getpid(),
                                                                       gsub(" |:|-","", Sys.time())),
                                              format='raster', overwrite = TRUE)
              model@binary <- writeRaster(model@binary, paste0(tmppath,
                                                               "/bin", model.name, Sys.getpid(),
                                                               gsub(" |:|-","", Sys.time())),
                                          format='raster', overwrite = TRUE)
            }
          }
          modelrep <- c(modelrep,model)
        }
        return(modelrep)
        }

        models <- unlist(models)

        if (verbose) {
          cat(paste(Sys.time(),"Finished modelling", paste0(algorithms,collapse="/"), "in parallel \n\n"))
        }
      } # end parmode algorithms

    stopCluster(cl)
    if(verbose){
      cat("Closed clusters")
    }

  } else {
  for (i in seq_len(length(algorithms))) {
    for (j in 1:rep) {
      model.name <- paste0(algorithms[i], ".", sprintf(paste0("%0",nchar(rep),"d"),j))
      if (verbose) {
        cat("Modelling :", model.name, "\n\n")
      }
      model <- try(modelling(algorithms[i], Occurrences, Env, Xcol = Xcol,
                             Ycol = Ycol, Pcol = Pcol, name = NULL, PA = PA,
                             cv = cv, cv.param = cv.param, final.fit.data = final.fit.data,
                             bin.thresh = bin.thresh, metric = metric, thresh = thresh,
                             axes.metric = axes.metric,
                             select = FALSE, select.metric = ensemble.metric,
                             select.thresh = ensemble.thresh,
                             verbose = verbose, GUI = GUI, ...))
      if (GUI) {
        incProgress(1/(length(algorithms) + 1), detail = paste(algorithms[i],
                                                               "SDM built"))
      }
      if (inherits(model, "try-error")) {
        if (verbose) {
          cat(model)
        }
      } else {
        ### not working yet
        if (!isFALSE(tmp) & !is.null(model)) {
          model@projection <- writeRaster(model@projection, paste0(tmppath,
                                                                   "/proba",model.name),
                                          format='raster', overwrite = TRUE)
          model@binary <- writeRaster(model@binary, paste0(tmppath,
                                                           "/bin",model.name),
                                      format='raster', overwrite = TRUE)
        }
        suppressWarnings({
          models[model.name] <- model
        })
      }
      if (verbose) {
        cat(sprintf("%s done for %s %s \n\n", model.name, spname, format(Sys.time(),"%Y-%m-%d %T")))
      }
    } # j replicates
  } # i algos
}
  # Ensemble modelling
  if (verbose) {
    cat(sprintf("#### Ensemble modelling with algorithms models for %s ##### %s \n\n",
                spname, format(Sys.time(), "%Y-%m-%d %T")))
  }
  algo <- list()
  for (i in seq_len(length(models))) {
    algo[[i]] <- models[[i]]
  }
  if (!is.null(name)) {
    algo["name"] <- name
  }
  algo["thresh"] <- thresh
  algo["uncertainty"] <- uncertainty
  algo["SDM.projections"] <- SDM.projections
  algo[["ensemble.metric"]] <- ensemble.metric
  algo[["ensemble.thresh"]] <- ensemble.thresh
  algo["weight"] <- weight
  algo["cores"] <- cores
  algo["verbose"] <- verbose
  esdm <- do.call(ensemble, algo)
  if (verbose) {
    cat(sprintf("Ensemble modelling done for %s %s \n\n", spname, format(Sys.time(),
                                                                         "%Y-%m-%d %T")))
  }
  if (GUI) {
    incProgress(1/(length(algorithms) + 1), detail = "Ensemble SDM built")
  }

  if (!is.null(esdm)) {
    # Parameters
    text.algorithms <- character()
    for (i in seq_len(length(algorithms))) {
      text.algorithms <- paste0(text.algorithms, ".", algorithms[i])
    }
    esdm@parameters$algorithms <- text.algorithms
    esdm@parameters$rep <- rep

    # Saving
    if (save) {
      if (verbose) {
        cat("#### Saving ##### \n\n")
      }
      save.esdm(esdm, path = path)
    }
  }

  # Removing tmp
  if (!isFALSE(tmp)) {
    unlink(paste0(tmppath), recursive = TRUE, force = TRUE)
  }
  rm(list = ls()[-which(ls() == "esdm")])
  gc()
  return(esdm)
}
