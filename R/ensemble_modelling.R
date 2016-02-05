#' @include modelling.R ensemble.R Ensemble.SDM.R checkargs.R
#' @importFrom shiny incProgress
#' @importFrom raster writeRaster
NULL

#'Build an ensemble SDM that assembles multiple algorithms
#'
#'This is a function to build an ensemble SDM that assembles multiple algorithms
#'for a single species. The function takes as inputs an occurrence data frame
#'made of presence/absence or presence-only records and a raster object for data
#'extraction and projection. The function returns an S4
#'\linkS4class{Ensemble.SDM} class object containing the habitat suitability
#'map, the binary map, the between-algorithm variance map and the associated
#'evaluation tables (model evaluation, algorithm evaluation, algorithm
#'correlation matrix and variable importance).
#'
#'@param algorithms character. Choice of the algorithm(s) to be run (see details
#'  below).
#'@param Occurrences data frame. Occurrence table (can be processed first by
#'  \code{\link{load_occ}}).
#'@param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurrence table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence, by setting presence to 1 and
#'  absence to 0. If NULL presence-only dataset is assumed.
#'@param rep integer. Number of repetitions for each algorithm.
#'@param name character. Optional name given to the final Ensemble.SDM produced
#'  (by default 'Ensemble.SDM').
#'@param save logical. If set to true, the ensemble SDM is automatically saved.
#'@param path character. If save is true, the path to the directory in which the
#'  ensemble SDM will be saved.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only dataset. If PA is NULL, recommended PA selection
#'  strategy is used depending on the algorithm (see details below).
#'@param cv character. Method of cross-validation used to evaluate the ensemble
#'  SDM (see details below).
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the ensemble SDM (see details below).
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate variable relative importance (see
#'  details below).
#'@param uncertainty logical. If set to true, generates an uncertainty map and
#'  an algorithm correlation matrix.
#'@param tmp logical. If set to true, the habitat suitability map of each
#'  algorithm is saved in a temporary file to release memory. But beware: if you
#'  close R, temporary files will be destroyed. To avoid any loss you can save
#'  your ensemble SDM with \code{\link{save.model}}. Depending on number,
#'  resolution and extent of models, temporary files can take a lot of disk
#'  space. Temporary files are written in R environment temporary folder.
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'@param weight logical. Choose whether or not you want the SDMs to be weighted
#'  using the selection metric or, alternatively, the mean of the selection
#'  metrics.
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'@param ... additional parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Ensemble.SDM} class object viewable with the
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
#'  plot (receiving operating characteristic curve) and the upper left corner
#'  (1,1).}\item{axes.metric}{Choice of the metric used to evaluate the variable
#'  relative importance (difference between a full model and one with each
#'  variable successively omitted): \strong{Pearson} (computes a simple
#'  Pearson's correlation \emph{r} between predictions of the full model and the
#'  one without a variable, and returns the score \emph{1-r}: the highest the
#'  value, the more influence the variable has on the model), \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{ensemble.metric}{Ensemble metric(s) sed to select SDMs: \strong{AUC},
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \item{"..."}{See algorithm in detail section} }
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
#'  be more or less time- and memory-consuming.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occurrences = subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#'
#' # ensemble SDM building
#' ESDM = ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 1,
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
#'@export
ensemble_modelling = function(algorithms,
                              # Modelling data input
                              Occurrences, Env,
                              # Occurrences reading
                              Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                              # Model creation
                              rep = 10, name = NULL, save = F, path = getwd(),
                              # Pseudo-absences definition
                              PA = NULL,
                              # Evaluation parameters
                              cv = 'holdout', cv.param = c(0.7,1), thresh = 1001, metric = 'SES',
                              axes.metric = 'Pearson', uncertainty = T, tmp = F,
                              # Assembling parameters
                              ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = T,
                              # Informations parameters
                              verbose = T, GUI = F,
                              # Modelling parameters
                              ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, rep = rep, name = name, save = save,
             path = path, PA = PA, cv = cv, cv.param = cv.param, thresh = thresh,
             metric = metric, axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
             ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, verbose = verbose, GUI = GUI)

  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM')
  if ('all' %in% algorithms) {algorithms = available.algo}
  for (i in 1:length(algorithms)) {
    if(!(algorithms[[i]] %in% available.algo)) {stop(algorithms[[i]],' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, SVM')}}
  if (tmp) {
    tmppath = get("tmpdir",envir=.PkgEnv)
    if (!("/.models" %in% list.dirs(path))) (dir.create(paste0(tmppath,'/.models')))
  }

  # Algorithms models creation
  if(is.null(name)){
    spname = 'species'
  } else {
    spname = name
  }
  if(verbose){cat(sprintf('#### Algorithms models creation for %s ##### %s \n\n', spname, format(Sys.time(), '%Y-%m-%d %T')))}
  models = list()
  for (i in 1:length(algorithms)) {
    for (j in 1:rep) {
      model.name = paste0(algorithms[i],'.',j)
      if(verbose){cat('Modelling :', model.name, '\n\n')}
      model = try(modelling(algorithms[i], Occurrences, Env, Xcol = Xcol, Ycol = Ycol, Pcol = Pcol,
                            name = NULL, PA = PA, cv = cv, cv.param = cv.param, thresh = thresh,
                            metric = metric, axes.metric =axes.metric, select = F,
                            select.metric = ensemble.metric, select.thresh = ensemble.thresh,
                            verbose = verbose, GUI = GUI, ...))
      if(GUI) {incProgress(1/(length(algorithms)+1),
                           detail = paste(algorithms[i],'SDM built'))}
      if (inherits(model, "try-error")) {if(verbose){cat(model)}} else {
        if (tmp) {model@projection = writeRaster(model@projection, paste0(tmppath,'/.models/',j,model.name), overwrite = T)}
        models[model.name] = model
      }
      if(verbose){cat(sprintf('%s done for %s %s \n\n', model.name, spname, format(Sys.time(), '%Y-%m-%d %T')))}
    }
  }

  # Ensemble modelling
  if(verbose){cat(sprintf('#### Ensemble modelling with algorithms models for %s ##### %s \n\n', spname, format(Sys.time(), '%Y-%m-%d %T')))}
  algo = list()
  for (i in 1:length(models)) {algo[[i]] = models[[i]]}
  if (!is.null(name)) {algo['name'] = name}
  algo['thresh'] = thresh
  algo['uncertainty'] = uncertainty
  algo[['ensemble.metric']] = ensemble.metric
  algo[['ensemble.thresh']] = ensemble.thresh
  algo['weight'] = weight
  enm = do.call(ensemble, algo)
  if(verbose){cat(sprintf('Ensemble modelling done for %s %s \n\n', spname, format(Sys.time(), '%Y-%m-%d %T')))}
  if(GUI) {incProgress(1/(length(algorithms)+1), detail = 'Ensemble SDM built')}

  if(!is.null(enm)) {
    # Parameters
    text.algorithms = character()
    for (i in 1:length(algorithms)) {text.algorithms = paste0(text.algorithms,'.',algorithms[i])}
    enm@parameters$algorithms = text.algorithms
    enm@parameters$rep = rep

    # Saving
    if(save) {
      if(verbose){cat('#### Saving ##### \n\n')}
      save.enm(enm, path = path)
    }
  }

  # Removing tmp
  if (tmp) {unlink('./.models', recursive = T, force = T)}
  rm(list = ls()[-which(ls() == 'enm')])
  gc()
  return(enm)
}
