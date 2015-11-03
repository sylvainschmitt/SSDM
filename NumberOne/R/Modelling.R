#' @include Algorithm.Niche.Model.R evaluate.R
#' @importFrom raster stack writeRaster
NULL

#'One algorithm specie distribution modelling
#'
#'This is a function for modelling one specie distribution with one algorithm.
#'It takes in inputs an occurences data frame made with presence/absence or
#'presence-only records and a raster objects for data extractions and
#'projections. It returns an S4 \linkS4class{Algorithm.Niche.Model} class object
#'containing the habitat suitability map and all evaluation tables comming with
#'(model evaluation and variables importance).
#'
#'@param algorithm character. Choice of the algorithm for the modelling (see
#'  details below).
#'@param Occurences data frame. Occurences table (can be treated first by
#'  load.occ).
#'@param Env raster object. Environnment raster object (can be treated first by
#'  load.env).
#'@param Xcol character. Name of the occurences table column containing Latitude
#'  or X coordinates.
#'@param Ycol character. Name of the occurences table column containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the occurences table column containing presence
#'  or absence value, if NULL presence-only data set is assumed.
#'@param name character. Optionnal name given to the final Algorithm.Niche.Model
#'  producted.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence only occurences, if PA is NULL recommended PA is used
#'  depending on the algorithms (see details below).
#'@param train.frac float [0,1]. training fraction used to train the model
#'  (1-train.frac being the evaluating fraction of the data set used to evaluate
#'  the model), in case of data split evaluation method (see details below).
#'@param thresh numeric. binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'@param metric character. Metric used to compute the binary map threshold (see
#'  details below.)
#'@param axes.metric Metric used to evaluate the variables relative importance
#'  in percent (see details below).
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Algorithm.Niche.Model} Class object viewable with
#'  \code{\link{plot}} method
#'
#'@details \describe{ \item{algorithm}{Currently available algorithms are
#'  Generalized linear models (\strong{GLM}), Generalized additive models
#'  (\strong{GAM}), Multivaraite adaptative regression splines (\strong{MARS}),
#'  Generalized boosted regressions models (\strong{GBM}), Classification tree
#'  analysis (\strong{CTA}), Random forests (\strong{RF}), Maximum Entropy
#'  (\strong{MAXENT}), Artificial Neural Network (\strong{ANN}, and Support
#'  vector machines (\strong{SVM}). Each algorithm have his own parameters
#'  settable with the (\strong{...}, see each algorithm section below to set
#'  theim.)} \item{PA}{list with two values : \strong{nb} number of pseudo
#'  absence selected, and \strong{strat} strategy used for pseudo-absence
#'  selection : either random selection either disk selection. We set to default
#'  the Barbet and Massin 200X recommendation (link).} \item{metric}{Choice of
#'  the metric used to compute binary map threshold and confusion matrix :
#'  \strong{Kappa} maximizes the model Kappa value, \strong{CCR} maximizes the
#'  correct predicted observations proportion, \strong{TSS} (True Skill
#'  Statistic) maximizes the sensitivity and specificity sum, \strong{SES} using
#'  the sensitivty specificity equality, \strong{LW} using the lowest occurence
#'  prediction probability, \strong{ROC} minimizing the distance between the ROC
#'  plot (receiving operative curve) and the upper left coin (1,1).}
#'  \item{axes.metric}{Choice of the metric used to evaluate the variables
#'  relative importance in percent (variation of the model evaluation without
#'  this axis) : \strong{Pearson} pearson correlation coefficient, \strong{AUC}
#'  area under the receiving operating curve (ROC), \strong{Kappa},
#'  \strong{omission.rate} omission rate, \strong{sensitivity},
#'  \strong{specificity}, and \strong{prop.correct} correct predicted
#'  occurences proportion.} \item{...}{See algorithm in detail section} }
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
#'  be more or less time and memory consuming. You can use the function estimate
#'  to have an estimation of both time and memory needed in your case with your
#'  hardware.
#'
#' @examples
#'\dontrun{
#' Modelling('GLM', Occurences, Env)
#'}
#'
#'@seealso \code{\link{Ensemble.Modelling}} for specie distribution ensemble
#'  modelling with multiple algorithms, \code{\link{Stack.Modelling}} for stack
#'  species ensemble modelling with multiple algorithms and species
#'
#'@export
Modelling = function(algorithm,
                     # Modelling data input
                     Occurences, Env,
                     # Occurences reading
                     Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL,
                     # Model creation
                     name = NULL,
                     # Pseudo-absences definition
                     PA = NULL, train.frac = 0.7,
                     # Evaluation parameters
                     thresh = 1001, metric = 'SES', axes.metric = 'AUC',
                     # Modelling parameters
                     ...) {
  # Test if algorithm is available
  available.algo = c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM')
  if (algorithm == 'all') {algorithm = available.algo}
  if(!(algorithm %in% available.algo)) {stop(algorithm,' is still not available, please use one of those : GLM, GAM, MARS, GBM, CTA, RF, MAXENT, ANN, SVM')}

  # Empty Algorithm niche model object creation
  model = Algorithm.Niche.Model(algorithm)
  if (!is.null(name)) {name = paste0(name,'.')}
  model@name = paste0(name,algorithm,'.Niche.Model')
  model@parameters$data = 'presence/absence data set'
  model@parameters$metric = metric

  cat('Data check ... \n')
  # Occurences data input test | Data frame needed
  if (is.matrix(Occurences)) {Occurences = data.frame(occurences)}
  if (!is.data.frame(Occurences)) {stop('Occurences data set is not a data frame or a matrix')}
  if ((Xcol %in% names(Occurences)) == F) {stop('X column is not well defined')}
  if ((Ycol %in% names(Occurences)) == F) {stop('Y column is not well defined')}
  if (is.null(Pcol)) {
    PO = T # Presence only
    cat('No presence column, presence-only data set is supposed.\n')
    model@parameters$data = 'presence-only data set'
  } else if ((Pcol %in% names(Occurences)) == F) {stop('Presence column is not well defined')}
  if (!is.null(PA)) {PO = T}
  if (PO) {cat('Pseudo-absence selection will be computed.\n')}
  data = data.frame(X = Occurences[which(names(Occurences) == Xcol)], Y = Occurences[which(names(Occurences) == Ycol)])
  names(data) = c('X','Y')
  if (PO) {data$Presence = 1} else {data$Presence = Occurences[which(names(Occurences == Pcol))]}
  data$Train = F
  data$Train[sample.int(length(data$Presence), round(length(data$Presence)*train.frac))] = T

  # Environment data input test | RasterStack needed
  if (is.raster(Env)) {Env = stack(Env)}
  if (!inherits(Env, 'RasterStack')) {stop('Environment data set is not a raster or a raster stack')}
  cat('   done. \n\n')

  # Pseudo - absences selection
  cat('Pseudo absence selection... \n')
  model@data = data
  if (PO) {
    model = PA.select(model, Env, PA, train.frac)
    model@parameters['PA'] = T}
  model = data.values(model, Env)
  cat('   done. \n\n')

  # Projection
  cat('Model projection...')
  model = project(model, Env, ...)
  cat('   done. \n\n')

  # Evaluation
  cat('Model evaluation...\n')
  model = evaluate(model, thresh, metric)
  cat('   done. \n\n')

  # Evaluation
  cat('Model axes contribution evaluation...\n')
  model = evaluate.axes(model, thresh, Env, axes.metric, ...)
  cat('   done. \n\n')

  return(model)
}
