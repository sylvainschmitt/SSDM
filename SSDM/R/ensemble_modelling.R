#' @include modelling.R ensemble.R Ensemble.SDM.R checkargs.R
#' @importFrom raster writeRaster
NULL

#'Build an an ensemble SDM that assembles multiple algorithms
#'
#'This is a function to build an ensemble SDM that assembles multiple algorithms
#'for a single species. The function takes as inputs an occurrence data frame
#'made of presence/absence or presence-only records and a raster object for data
#'extraction and projection. The function returns an S4
#'\linkS4class{Ensemble.SDM} class object containing the habitat suitability
#'map, the binary map, and the uncertainty map (based on the between-algorithm
#'variance of habitat suitability maps) and the associated evaluation tables
#'(model evaluation, algorithm evaluation, algorithm correlation matrix and
#'variable importance).
#'
#'@param algorithms character. Choice of the algorithm(s) to be run (see details
#'  below).
#'@param Occurrences data frame. Occurrence table (can be processed first by
#'  \code{\link{load_occ}}).
#'@param Env raster object. Raster object of environmental variable (can be
#'  processed first by \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurences table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  wether a line is a presence or an absence. If NULL presence-only data set is
#'  assumed.
#'@param rep integer. Number of repetitions for each algorithm.
#'@param name character. Optional name given to the final Ensemble.SDM producted
#'  (by default 'Ensemble.SDM').
#'@param save logical. If true the ensemble SDM is automatically saved.
#'@param directory character. If save is true, the name of the directory in
#'  which the ensemble SDM will be saved.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only data. If PA is NULL recommended PA is used
#'  depending on the algorithm (see details below).
#'@param cv character. Method of cross-validation used to evaluate the ensemble
#'  SDM (see details below).
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the ensemble SDM (see details below).
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 & 1. The higher it is the more accurate
#'  is the threshold but the longer is the modelling evaluation step (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate variable relative importance (see
#'  details below).
#'@param uncertainty logical. If true generates an uncertainty map and algorithm
#'  correlation matrix are computed.
#'@param tmp logical. If true the habitat suitability map of each algorithm is
#'  saved in a temporary file to release memory. But beware: if you close R,
#'  temporary files will be destroyed. To avoid any loss you can save your
#'  ensemble SDM with \code{\link{save.model}}.
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'@param weight logical. Choose wether or not you want the SDMs to be weighted
#'  using the slection metric or, alternatively, the mean of the selection
#'  metrics.
#'@param verbose logical. If true allow the function to print text in the
#'  console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Ensemble.SDM} class object viewable with
#'  \code{\link{plot.model}} function.
#'
#'@details \describe{ \item{algorithms}{'all' allows you to call directly all
#'  available algorithms. Currently, available algorithms include Generalized
#'  linear model (\strong{GLM}), Generalized additive model (\strong{GAM}),
#'  Multivaraite adaptative regression splines (\strong{MARS}), Generalized
#'  boosted regressions model (\strong{GBM}), Classification tree analysis
#'  (\strong{CTA}), Random forest (\strong{RF}), Maximum Entropy
#'  (\strong{MAXENT}), Artificial Neural Network (\strong{ANN}), and Support
#'  vector machines (\strong{SVM}). Each algorithm has its own parameters
#'  settable with the (\strong{...}, see each algorithm section below to set
#'  their parameters.)} \item{"PA"}{list with two values : \strong{nb} number of
#'  pseudo-absence selected, and \strong{strat} strategy used for pseudo-absence
#'  selection: either random selection or disk selection. We set to default
#'  recommendation from Barbet-Massin et al. (2012) (see references).}
#'  \item{cv}{\strong{Cross validation} method used to split the occurence
#'  dataset for evaluation: \strong{holdout} data are partitioned in a training
#'  set and an evaluating set using a fraction (\emph{cv.param[1]}) and the
#'  operation can be repeated (\emph{cv.param[2]}), \strong{k-folds} data are
#'  partitioned in k folds successively being the evaluating set regarding a
#'  parameter k (\emph{cv.param[1]}) and the operation can be repeated
#'  (\emph{cv.param[2]}), \strong{LOO} (Leave One Out) each point are
#'  successively take as evaluation data.} \item{metric}{Choice of the metric
#'  used to compute the binary map threshold and the confusion matrix (by
#'  default SES as recommended by Liu et al. (2005), see references below):
#'  \strong{Kappa} maximizes the Kappa, \strong{CCR} maximizes the proportion of
#'  correctly predicted observations, \strong{TSS} (True Skill Statistic)
#'  maximizes the sum of sensitivity and specificity, \strong{SES} using the
#'  sensitivty-specificity equality, \strong{LW} using the lowest occurrence
#'  prediction probability, \strong{ROC} minimizing the distance between the ROC
#'  plot (receiving operating curve) and the upper left corner
#'  (1,1).}\item{axes.metric}{Choice of the metric used to evaluate the variable
#'  relative importance in percent (variation of the model evaluation without
#'  this axis): \strong{Pearson} Pearson's correlation coefficient, \strong{AUC}
#'  area under the receiving operating characteristic (ROC) curve,
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} proportion of correctly predicted occurrences.}
#'  \item{ensemble.metric}{Ensemble metric (metric used to compute the SDMs
#'  selection): \strong{AUC} area under the receiving operating characteristic
#'  (ROC) curve, \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} proportion of correctly predicted occurrences.}
#'  \item{"..."}{See algorithm in detail section} }
#'
#'@section Generalized linear model (\strong{GLM}) : Uses the glm function from
#'  the package 'stats', you can set the following parameters (see
#'  \code{\link[stats]{glm}} for more details): \describe{
#'  \item{test}{character. Test used to evaluate the model, default 'AIC'.}
#'  \item{epsilon}{numeric. Epsilon value used to fit the model, default
#'  10e-08.} \item{maxit}{numeric. Maximum number of iterations allowed to fit
#'  the SDM, default 500.} }
#'
#'@section Generalized additive model (\strong{GAM}) : Uses the gam function
#'  from the package 'mgcv', you can set the following parameters (see
#'  \code{\link[mgcv]{gam}} for more details): \describe{ \item{test}{character.
#'  Test used to evaluate the model, default 'AIC'.} \item{epsilon}{numeric.
#'  Epsilon value used to fit the model, default 10e-08.} \item{maxit}{numeric.
#'  Maximum number of iterations allowed to fit the SDM, default 500.} }
#'
#'@section Multivaraite adaptative regression splines (\strong{MARS}) : Uses the
#'  earth function from the package 'earth', you can set the following
#'  parameters (see \code{\link[earth]{earth}} for more details): \describe{
#'  \item{degree}{integer. Number of interaction degrees allowed in the SDM,
#'  default 2.} }
#'
#'@section Generalized boosted regressions model (\strong{GBM}) : Uses the gbm
#'  function from the package 'gbm,' you can set the following parameters (see
#'  \code{\link[gbm]{gbm}} for more details): \describe{ \item{trees}{integer.
#'  Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leave of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} \item{thresh.shrink}{integer. Tree shrinkage coefficient,
#'  default 1e-03.} }
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the rpart function
#'  from the package 'rpart', you can set the following parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leave of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} }
#'
#'@section Random Forest (\strong{RF}) : Uses the randomForest function from the
#'  package 'randomForest', you can set those parameters (see
#'  \code{\link[randomForest]{randomForest}} for more details): \describe{
#'  \item{trees}{integer. Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leave of trees, default 1.} }
#'
#'@section Maximum Entropy (\strong{MAXENT}) : It uses the maxent function from
#'  the package 'dismo'. Take care to have correctly install the maxent.jar file
#'  in the folder dismo in your library folder (see \code{\link[dismo]{maxent}}
#'  for more details).
#'
#'@section Artificial Neural Network (\strong{ANN}) : Uses the nnet function
#'  from the package 'nnet', you can set the following parameters (see
#'  \code{\link[nnet]{nnet}} for more details): \describe{ \item{maxit}{integer.
#'  Maximum number of iterations, default 500.} }
#'
#'@section Support vector machines (\strong{SVM}) : Uses the svm function from
#'  the package 'e1071', you can set the following parameters (see
#'  \code{\link[e1071]{svm}} for more details): \describe{ \item{epsilon}{float.
#'  Epsilon parameter in the insensitive loss function, default 1e-08.}
#'  \item{cv}{integer. Number of cross-validations, default 3.} }
#'
#'@section Warning : Depending on the raster object resolution the computing can
#'  be more or less time- and memory-consuming.
#'
#' @examples
#'\dontrun{
#' ensemble.modelling('all', Occurrences, Env)
#'}
#'
#'@seealso \code{\link{modelling}} to build SDMs with a single algorithm,
#'  \code{\link{stack_modelling}} tu build SSDMs.
#'
#'@references M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012)
#'  "Selecting pseudo-absences for species distribution models: how, where and
#'  how many?" \emph{Methods Ecology and Evolution} 3(2):327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting
#'  thresholds of occurrence in the prediction of species distributions."
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'@export
ensemble_modelling = function(algorithms,
                              # Modelling data input
                              Occurrences, Env,
                              # Occurrences reading
                              Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                              # Model creation
                              rep = 1, name = NULL, save = F, directory = getwd(),
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
             directory = directory, PA = PA, cv = cv, cv.param = cv.param, thresh = thresh,
             metric = metric, axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
             ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, verbose = verbose, GUI = GUI)

  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM')
  if (algorithms == 'all') {algorithms = available.algo}
  for (i in 1:length(algorithms)) {
    if(!(algorithms[[i]] %in% available.algo)) {stop(algorithms[[i]],' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, SVM')}}
  if (tmp) {
    path = get("tmpdir",envir=.PkgEnv)
    if (!("/.models" %in% list.dirs(path))) (dir.create(paste0(path,'/.models')))
  }

  # Algorithms models creation
  if(verbose){cat('#### Algorithms models creation ##### \n\n')}
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
        if (tmp) {model@projection = writeRaster(model@projection, paste0(path,'/.models/',j,model.name), overwrite = T)}
        models[model.name] = model
      }
      if(verbose){cat('\n\n')}
    }
  }

  # Ensemble modelling
  if(verbose){cat('#### Ensemble modelling with algorithms models ##### \n\n')}
  algo = list()
  for (i in 1:length(models)) {algo[[i]] = models[[i]]}
  if (!is.null(name)) {algo['name'] = name}
  algo['thresh'] = thresh
  algo['uncertainty'] = uncertainty
  algo[['ensemble.metric']] = ensemble.metric
  algo[['ensemble.thresh']] = ensemble.thresh
  algo['weight'] = weight
  enm = do.call(ensemble, algo)
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
      save.enm(enm, directory = directory)
    }
  }

  # Removing tmp
  if (tmp) {unlink('./.models', recursive = T, force = T)}

  return(enm)
}
