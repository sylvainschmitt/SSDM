#' @include Algorithm.SDM.R checkargs.R
#' @importFrom raster stack writeRaster
NULL

#'Build an SDM using one algorithm
#'
#'This is a function for modelling one specie distribution with one algorithm.
#'It takes in inputs an occurences data frame made with presence/absence or
#'presence-only records and a raster objects for data extractions and
#'projections. It returns an S4 \linkS4class{Algorithm.SDM} class object
#'containing the habitat suitability map, and the binary map and the evaluation
#'table.
#'
#'@param algorithm character. Choice of the algorithm (see details below).
#'@param Occurrences data frame. Occurrences table (can be treated first by
#'  \code{\link{load_occ}}).
#'@param Env raster object. Raster object of environmental variable (can be
#'  treated first by \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurences table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurences table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurences table specifying
#'  wether a line is a presence or an absence. If NULL presence-only data set is
#'  assumed.
#'@param name character. Optional name given to the final SDM producted.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only data. If PA is NULL recommended PA is used
#'  depending on the algorithms (see details below).
#'@param cv character. Method of cross-validation used to evaluate the SDM (see
#'  details below).
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the SDM (see details below).
#'@param select logical. If true models are evaluated before being projected on
#'  a raster, and don't keep if they don't match selection criteria. (see
#'  details below).
#'@param select.metric character. Metric used to pre-select SDMs that reach a
#'  sufficient quality.
#'@param select.thresh numeric. Threshold associated with the metric used to
#'  compute the selection.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 & 1. The higher it is the more accurate
#'  is the threshold but the longer is the modelling evaluation step (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate variable relative importance (see
#'  details below).
#'@param verbose logical. If true allow the function to print text in the
#'  console
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface) !
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Algorithm.SDM} Class object viewable with
#'  \code{\link{plot.model}} method
#'
#'@details \describe{ \item{algorithm}{'all' allows you to call directly all
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
#'  \item{cv}{\strong{Cross validation} method used to split the occurences
#'  dataset for evaluation: \strong{holdout} data are partitioned in a training
#'  set and an evaluating set using a fraction (\emph{cv.param[1]}) and the
#'  operation can be repeated (\emph{cv.param[2]}), \strong{k-folds} data are
#'  partitioned in k folds successively being the evaluating set regarding a
#'  parameter k (\emph{cv.param[1]}) and the operation can be repeated
#'  (\emph{cv.param[2]}), \strong{LOO} (Leave One Out) each point are
#'  successively take as evaluation data.} \item{metric}{Choice of the metric
#'  used to compute the binary map threshold and the confusion matrix (by
#'  default SES as recommended by Liu et al. (2005),see references below):
#'  \strong{Kappa} maximizes the Kappa, \strong{CCR} maximizes the proportion of
#'  correctly predicted observations, \strong{TSS} (True Skill Statistic)
#'  maximizes the sum of sensitivity and specificity, \strong{SES} using the
#'  sensitivty-specificity equality, \strong{LW} using the lowest occurence
#'  prediction probability, \strong{ROC} minimizing the distance between the ROC
#'  plot (receiving operative curve) and the upper left corner
#'  (1,1).}\item{axes.metric}{Choice of the metric used to evaluate the variable
#'  relative importance in percent (variation of the model evaluation without
#'  this axis): \strong{Pearson} Pearson's correlation coefficient, \strong{AUC}
#'  area under the receiving operating characteristic (ROC) curve,
#'  \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} proportion of correctly predicted occurences.}
#'  \item{select.metric}{Ensemble metric (metric used to compute the SDMs
#'  selection): \strong{AUC} area under the receiving operating characteristic
#'  (ROC) curve, \strong{Kappa}, \strong{sensitivity}, \strong{specificity}, and
#'  \strong{prop.correct} proportion of correctly predicted occurences.}
#'  \item{"..."}{See algorithm in detail section} }
#'
#'@section Generalized linear models (\strong{GLM}) : Uses the glm function from
#'  the package 'stats', you can set the following parameters (see
#'  \code{\link[stats]{glm}} for more details): \describe{
#'  \item{test}{character. Test used to evaluate the model, default 'AIC'.}
#'  \item{epsilon}{numeric. Epsilon value used to fit the model, default
#'  10e-08.} \item{maxit}{numeric. Maximum number of iterations allowed to fit
#'  the SDM, default 500.} }
#'
#'@section Generalized additive models (\strong{GAM}) : Uses the gam function
#'  from the package 'mgcv', you can set the following parameters (see
#'  \code{\link[mgcv]{gam}} for more details): \describe{ \item{test}{character.
#'  Test used to evaluate the model, default 'AIC'.} \item{epsilon}{numeric.
#'  Epsilon value used to fit the model, default 10e-08.} \item{maxit}{numeric.
#'  Maximum number of iterations allowed to fit the SDM, default 500.} }
#'
#'@section Multivaraite adaptative regression splines (\strong{MARS}) : Uses the
#'  earth function from the package 'earth', you can set the following
#'  parameters (see \code{\link[earth]{earth}} for more details): \describe{
#'  \item{degree}{integer. Number of interactions degrees allowed in the SDM,
#'  default 2.} }
#'
#'@section Generalized boosted regressions models (\strong{GBM}) : Uses the gbm
#'  function from the package 'gbm,' you can set the following parameters (see
#'  \code{\link[gbm]{gbm}} for more details): \describe{ \item{trees}{integer.
#'  Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} \item{thresh.shrink}{integer. Tree shrinkage coefficient,
#'  default 1e-03.} }
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the rpart function
#'  from the package 'rpart', you can set the following parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} \item{cv}{integer. Number of cross-validations,
#'  default 3.} }
#'
#'@section Random Forest (\strong{RF}) : Uses the randomForest function from the
#'  package 'randomForest', you can set those parameters (see
#'  \code{\link[randomForest]{randomForest}} for more details): \describe{
#'  \item{trees}{integer. Number of trees used in the model, default 2500.}
#'  \item{final.leave}{integer. Minimum of observations allowed in the final
#'  leaves of trees, default 1.} }
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
#' modelling('GLM', Occurrences, Env)
#'}
#'
#'@seealso \code{\link{ensemble_modelling}} for ensemble SDMs with multiple
#'  algorithms, \code{\link{stack_modelling}} for SSDMs.
#'
#'@references M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012)
#'  "Selecting pseudo-absences for species distribution models: how, where and
#'  how many?" \emph{Methods Ecology and Evolution} 3(2):327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting thresholds
#'  of occurrence in the prediction of species distributions." \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'@export
modelling = function(algorithm,
                     # Modelling data input
                     Occurrences, Env,
                     # Occurrences reading
                     Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                     # Model creation
                     name = NULL,
                     # Pseudo-absences definition
                     PA = NULL,
                     # Evaluation parameters
                     cv = 'holdout', cv.param = c(0.7,2), thresh = 1001, metric = 'SES', axes.metric = 'Pearson',
                     # Selection parameters
                     select = F, select.metric = c('AUC'), select.thresh = c(0.75),
                     # Informations parameters
                     verbose = T, GUI = F,
                     # Modelling parameters
                     ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, name = name, PA = PA, cv = cv, cv.param = cv.param,
             thresh = thresh, metric = metric, axes.metric = axes.metric, select = select,
             select.metric = select.metric, select.thresh = select.thresh, verbose = verbose, GUI = GUI)

  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM')
  if (algorithm == 'all') {algorithm = available.algo}
  if(!(algorithm %in% available.algo)) {stop(algorithm,' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, SVM')}

  # Empty Algorithm niche model object creation
  model = Algorithm.SDM(algorithm)
  if (!is.null(name)) {name = paste0(name,'.')}
  model@name = paste0(name,algorithm,'.SDM')
  model@parameters$data = 'presence/absence data set'
  model@parameters$metric = metric

  cat('Data check ... \n')
  # Occurrences data input test | Data frame needed
  if (is.matrix(Occurrences)) {Occurrences = data.frame(Occurrences)}
  if (!is.data.frame(Occurrences)) {stop('Occurrences data set is not a data frame or a matrix')}
  if ((Xcol %in% names(Occurrences)) == F) {stop('X column is not well defined')}
  if ((Ycol %in% names(Occurrences)) == F) {stop('Y column is not well defined')}
  if (is.null(Pcol)) {
    PO = T # Presence only
    cat('No presence column, presence-only data set is supposed.\n')
    model@parameters$data = 'presence-only data set'
  } else if ((Pcol %in% names(Occurrences)) == F) {stop('Presence column is not well defined')}
  if (!is.null(PA)) {PO = T}
  if (PO) {cat('Pseudo-absence selection will be computed.\n')}
  data = data.frame(X = Occurrences[which(names(Occurrences) == Xcol)], Y = Occurrences[which(names(Occurrences) == Ycol)])
  names(data) = c('X','Y')
  if (PO) {data$Presence = 1} else {data$Presence = Occurrences[which(names(Occurrences == Pcol))]}

  # Environment data input test | RasterStack needed
  if (is.raster(Env)) {Env = stack(Env)}
  if (!inherits(Env, 'RasterStack')) {stop('Environment data set is not a raster or a raster stack')}
  cat('   done. \n\n')

  # Pseudo - absences selection
  cat('Pseudo absence selection... \n')
  model@data = data
  if (PO) {
    model = PA.select(model, Env, PA)
    model@parameters['PA'] = T}
  model = data.values(model, Env)
  cat('   done. \n\n')

  # Evaluation
  cat('Model evaluation...\n')
  model = evaluate(model, cv, cv.param, thresh, metric, Env, ...)
  cat('   done. \n\n')

  # Model selection
  test = T
  if(select) {
    for (j in 1:length(select.metric)) {
      if(model@evaluation[,which(names(model@evaluation) == select.metric[j])] < select.thresh[j]) {
        test = F
      }
    }
  }
  if(test){
    # Projection
    cat('Model projection...\n')
    model = project(model, Env, ...)
    cat('   done. \n\n')

    # Axes evaluation
    cat('Model axes contribution evaluation...\n')
    model = evaluate.axes(model, cv, cv.param, thresh, metric, axes.metric, Env, ...)
    cat('   done. \n\n')

    return(model)
  } else {
    cat('Model have been rejected, NULL is returned ! \n')
    return(NULL)
  }
}
