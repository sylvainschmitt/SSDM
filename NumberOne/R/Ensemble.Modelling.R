#' @include Modelling.R ensemble.R Ensemble.Niche.Model.R checkargs.R
#' @importFrom raster writeRaster
NULL

#'Multi-algorithms specie distribution ensemble modelling
#'
#'This is a function for modelling one specie distribution with several
#'algorithms and make an ensemble model of theim. It takes in inputs an
#'occurences data frame made with presence/absence or presence-only records and
#'a raster objects for data extractions and projections. It returns an S4
#'\linkS4class{Ensemble.Niche.Model} class object containing the habitat
#'suitability map, the binary map, and the uncertainty map based on the habitat
#'suitability map variance inter algorithms and all evaluation tables comming
#'with (model evaluation, algorithms evaluation, algorithms correlation matrix
#'and variables importance).
#'
#'@param algorithms character. Choice of the algorithm for the modelling (see
#'  details below).
#'@param Occurences data frame. Occurences table (can be treated first by
#'  \code{\link{load.occ}}).
#'@param Env raster object. Environnment raster object (can be treated first by
#'  \code{\link{load.var}}).
#'@param Xcol character. Name of the occurences table column containing Latitude
#'  or X coordinates.
#'@param Ycol character. Name of the occurences table column containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the occurences table column containing presence
#'  or absence value, if NULL presence-only data set is assumed.
#'@param rep integer. Number of repetition for each algorithm.
#'@param name character. Optionnal name given to the final Ensemble.Niche.Model
#'  producted.
#'@param save logical. If true the model is automatically saved.
#'@param directory character. If save is true, the name of the directory to save
#'  the model.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence only occurences, if PA is NULL recommended PA is used
#'  depending on the algorithms (see details below).
#'@param cv character. Type of cross-validation used to evaluate the model (see
#'  details below).
#'@param cv.param numeric. Parameters associated to the cross-validation used to
#'  evaluate the model (see details below).
#'@param thresh numeric. binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate the variables relative importance
#'  in percent (see details below).
#'@param uncertainity logical. If false uncertainity mapping and algorithms
#'  correlation matrix is not computed.
#'@param tmp logical. If true algorithms habitat suitability map is saved in
#'  temporary files to release memory. But beware, if you close R temporary
#'  files will be destructed. To avoid any loss you can save your model with
#'  \code{\link{save.model}}.
#'@param ensemble.metric character. Metric used to compute the selection among
#'  algorithms different models (see details below)
#'@param ensemble.thresh numeric. Threshold associated with the metric used to
#'  compute the selection.
#'@param weight logical. Choose if the model are weighted or not by the
#'  selection associated metrics mean.
#'@param verbose logical. If true allow the function to print text in the
#'  console
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface) !
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Ensemble.Niche.Model} Class object viewable with
#'  \code{\link{plot.model}} method
#'
#'@details \describe{ \item{algorithm}{'all' allows you to call directly all
#'  availables algorithms. Currently available algorithms are Generalized linear
#'  models (\strong{GLM}), Generalized additive models (\strong{GAM}),
#'  Multivaraite adaptative regression splines (\strong{MARS}), Generalized
#'  boosted regressions models (\strong{GBM}), Classification tree analysis
#'  (\strong{CTA}), Random forests (\strong{RF}), Maximum Entropy
#'  (\strong{MAXENT}), Artificial Neural Network (\strong{ANN}), and Support
#'  vector machines (\strong{SVM}). Each algorithm have his own parameters
#'  settable with the (\strong{...}, see each algorithm section below to set
#'  theim.)} \item{"PA"}{list with two values : \strong{nb} number of pseudo
#'  absence selected, and \strong{strat} strategy used for pseudo-absence
#'  selection : either random selection either disk selection. We set to default
#'  the Barbet-Massin 2012 recommendation (see references).}
#'  \item{cv}{\strong{Cross validation} method used for the evaluation among :
#'  \strong{holdout} data are partitionned in a training set and evaluating set
#'  regarding a fraction (\emph{cv.param[1]}) and the operation can be repeated
#'  (\emph{cv.param[2]}), \strong{k-fols} data are partitionned in k fold
#'  succesively being the evaluating set regarding a k parameter
#'  (\emph{cv.param[1]}) and the operation can be repeated (\emph{cv.param[2]}),
#'  \strong{LOO} (Leave One Out) each point are successively take as evaluation
#'  data.} \item{metric}{Choice of the metric used to compute binary map
#'  threshold and confusion matrix (by default SES as recommanded by Liu et al.
#'  2005,see references below): \strong{Kappa} maximizes the model Kappa value,
#'  \strong{CCR} maximizes the correct predicted observations proportion,
#'  \strong{TSS} (True Skill Statistic) maximizes the sensitivity and
#'  specificity sum, \strong{SES} using the sensitivty specificity equality,
#'  \strong{LW} using the lowest occurence prediction probability, \strong{ROC}
#'  minimizing the distance between the ROC plot (receiving operative curve) and
#'  the upper left coin (1,1).}\item{axes.metric}{Choice of the metric used to
#'  evaluate the variables relative importance in percent (variation of the
#'  model evaluation without this axis) : \strong{Pearson} pearson correlation
#'  coefficient, \strong{AUC} area under the receiving operating curve (ROC),
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} correct predicted occurences proportion.}
#'  \item{ensemble.metric}{Ensemble metric (Metric used to compute the selection
#'  among algorithms different models) can be choosed among those : \strong{AUC}
#'  area under the receiving operating curve (ROC), \strong{Kappa},
#'  \strong{sensitivity}, \strong{specificity}, and \strong{prop.correct}
#'  correct predicted occurences proportion.} \item{"..."}{See algorithm in
#'  detail section} }
#'
#'@section Generalized linear models (\strong{GLM}) : it uses the glm function
#'  from package stats, you can set those parameters (see
#'  \code{\link[stats]{glm}} for more details): \describe{
#'  \item{test}{character. Test used to evaluate the model, default 'AIC'.}
#'  \item{epsilon}{numeric. Epsilon value used to fit the model, default
#'  10e-08.} \item{maxit}{numeric. Maximimu iterations allowed to fit the model,
#'  default 500.} }
#'
#'@section Generalized additive models (\strong{GAM}) : It uses the gam function
#'  from package mgcv, you can set those parameters (see \code{\link[mgcv]{gam}}
#'  for more details): \describe{ \item{test}{character. Test used to evaluate
#'  the model, default 'AIC'.} \item{epsilon}{numeric. Epsilon value used to fit
#'  the model, default 10e-08.} \item{maxit}{numeric. Maximimu iterations
#'  allowed to fit the model, default 500.} }
#'
#'@section Multivaraite adaptative regression splines (\strong{MARS}) : It uses
#'  the earth function from package earth, you can set those parameters (see
#'  \code{\link[earth]{earth}} for more details): \describe{
#'  \item{degree}{integer. Number of interactions degrees allowed in the model,
#'  default 2.} }
#'
#'@section Generalized boosted regressions models (\strong{GBM}) : It uses the
#'  gbm function from package gbm, you can set those parameters (see
#'  \code{\link[gbm]{gbm}} for more details): \describe{ \item{trees}{integer.
#'  Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} \item{thresh.shrink}{integer. Trees shrinkage coefficient,
#'  default 1e-03.} }
#'
#'@section Classification tree analysis (\strong{CAT}) : It uses the rpart
#'  function from package rpart, you can set those parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} }
#'
#'@section Random Forest (\strong{RF}) : It uses the randomForest function from
#'  package randomForest, you can set those parameters (see
#'  \code{\link[randomForest]{randomForest}} for more details): \describe{
#'  \item{trees}{integer. Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} }
#'
#'@section Maximum Entropy (\strong{MAXENT}) : It uses the maxent function from
#'  package dismo. Take care to have correctly put the maxent.jar file in dismo
#'  folder in your library folder (see \code{\link[dismo]{maxent}} for more
#'  details).
#'
#'@section Artificial Neural Network (\strong{ANN}) : It uses the nnet function
#'  from package nnet, you can set those parameters (see
#'  \code{\link[nnet]{nnet}} for more details): \describe{ \item{maxit}{integer.
#'  Maximum number of iteration, default 500.} }
#'
#'@section Support vector machines (\strong{SVM}) : it uses the svm function
#'  from package e1071, you can set those parameters (see
#'  \code{\link[e1071]{svm}} for more details): \describe{ \item{epsilon}{float.
#'  Epsilon parameter in the insensitive loss function , default 1e-08.}
#'  \item{cv}{integer. Number of cross-validations, default 3.} }
#'
#'@section Warning : Depending on the raster object resolution the computing can
#'  be more or less time and memory consuming.
#'
#' @examples
#'\dontrun{
#' Ensemble.Modelling('all', Occurences, Env)
#'}
#'
#'@seealso \code{\link{Modelling}} for specie distribution modelling with one
#'  algorithm, \code{\link{Stack.Modelling}} for stack species ensemble
#'  modelling with multiple algorithms and species
#'
#'@references Barbet-Massin, M., Jiguet, F., Albert, C. H. & Thuiller, W.
#'  Selecting pseudo-absences for species distribution models: how, where and
#'  how many? Methods Ecol. Evol. 3, 327-338 (2012)
#'
#'  Liu,  C.  et  al.  2005.  Selecting  thresholds  of  occurrence  in  the
#'  prediction  of  species  distributions./ Ecography  28:  385 / 393.
#'
#'@export
Ensemble.Modelling = function(algorithms,
                              # Modelling data input
                              Occurences, Env,
                              # Occurences reading
                              Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                              # Model creation
                              rep = 1, name = NULL, save = F, directory = getwd(),
                              # Pseudo-absences definition
                              PA = NULL,
                              # Evaluation parameters
                              cv = 'holdout', cv.param = c(0.7,1), thresh = 1001, metric = 'SES',
                              axes.metric = 'Pearson', uncertainity = T, tmp = F,
                              # Assembling parameters
                              ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = T,
                              # Informations parameters
                              verbose = T, GUI = F,
                              # Modelling parameters
                              ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, rep = rep, name = name, save = save,
             directory = directory, PA = PA, cv = cv, cv.param = cv.param, thresh = thresh,
             metric = metric, axes.metric = axes.metric, uncertainity = uncertainity, tmp = tmp,
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
  cat('#### Algorithms models creation ##### \n\n')
  models = list()
  for (i in 1:length(algorithms)) {
    for (j in 1:rep) {
      model.name = paste0(algorithms[i],'.',j)
      cat('Modelling :', model.name, '\n\n')
      model = try(Modelling(algorithms[i], Occurences, Env, Xcol = Xcol, Ycol = Ycol, Pcol = Pcol,
                            name = NULL, PA = PA, cv = cv, cv.param = cv.param, thresh = thresh,
                            metric = metric, axes.metric =axes.metric, select = F,
                            select.metric = ensemble.metric, select.thresh = ensemble.thresh,...))
      if (inherits(model, "try-error")) {cat(model)} else {
        if (tmp) {model@projection = writeRaster(model@projection, paste0(path,'/.models/',j,model.name), overwrite = T)}
        models[model.name] = model
      }
      cat('\n\n')
    }
  }

  # Ensemble modelling
  cat('#### Ensemble modelling with algorithms models ##### \n\n')
  algo = list()
  for (i in 1:length(models)) {algo[[i]] = models[[i]]}
  if (!is.null(name)) {algo['name'] = name}
  algo['thresh'] = thresh
  algo['uncertainity'] = uncertainity
  algo[['ensemble.metric']] = ensemble.metric
  algo[['ensemble.thresh']] = ensemble.thresh
  algo['weight'] = weight
  enm = do.call(ensemble, algo)

  if(!is.null(enm)) {
    # Parameters
    text.algorithms = character()
    for (i in 1:length(algorithms)) {text.algorithms = paste0(text.algorithms,'.',algorithms[i])}
    enm@parameters$algorithms = text.algorithms
    enm@parameters$rep = rep

    # Saving
    if(save) {
      cat('#### Saving ##### \n\n')
      save.enm(enm, directory = directory)
    }
  }

  # Removing tmp
  if (tmp) {unlink('./.models', recursive = T, force = T)}

  return(enm)
}
