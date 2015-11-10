#' NumberOne: Stack species ensemble distribution modelling
#'
#' Number one is a package to realize stack species ensemble distribution
#' modelling. Beginning from simple specie distribution model with one
#' algorithm, you can then create an ensemble model with multiple algorithm
#' models from different algorithms and with repetition. Finally you can stack
#' all your species ensemble distribution models in a stack species model.
#' Models contains habitat suitability map, habitat binary map, uncertainty map,
#' evaluation, variables importance, algorithms evaluation and correlation, and
#' diversity map features.
#'
#' NumberOne package contain also a Global User Interface for easier utilisation
#' (\code{\link{NumberOneGUI}}).
#'
#' Number one provides 4 categories of functions (that you can find in details
#' below): Data preparation, Modelling main functions, Model main methods, Model
#' classes, and Miscellaneous.
#'
#' @section Data preparation: \describe{\item{\code{\link{load.occ}}}{Load
#'   occurences data set} \item{\code{\link{load.var}}}{Load environmental
#'   varaibles rasters} }
#'
#' @section Modelling main  functions: \describe{
#'   \item{\code{\link{Modelling}}}{To do one specie ditribution model with one
#'   algorithm} \item{\code{\link{Ensemble.Modelling}}}{To do one specie
#'   ensemble distribution model with multiple algorithms and repetitions}
#'   \item{\code{\link{Stack.Modelling}}}{To do stack species ensemble
#'   distribution model with multiple species, algorithms and repetitions}}
#'
#' @section Model main methods: \describe{
#'   \item{\code{\link{ensemble,Algorithm.Niche.Model-method}}}{Make an ensemble
#'   distribution model from multiple simple models}
#'   \item{\code{\link{stacking,Ensemble.Niche.Model-method}}}{Make a stack
#'   species ensemble distribution model from multiple ensemble distribution
#'   model}
#'   \item{\code{\link{update,Stack.Species.Ensemble.Niche.Model-method}}}{Update
#'    a previous stack with new occurences data}}
#'
#'
#' @section Model classes: \describe{
#'   \item{\code{\linkS4class{Algorithm.Niche.Model}}}{S4 class to represent a
#'   specie distribution model of one algorithm}
#'   \item{\code{\linkS4class{Ensemble.Niche.Model}}}{S4 class to represent a
#'   specie ensemble distribution model}
#'   \item{\code{\linkS4class{Stack.Species.Ensemble.Niche.Model}}}{S4 class to
#'   represent a stack species ensemble distribution model}}
#'
#' @section Miscellaneous: \describe{ \item{\code{\link{NumberOneGUI}}}{A global
#'   user interface for easier utilisation of the NumberOne package}
#'   \item{\code{\link{plot.model}}}{Plot models}
#'   \item{\code{\link{save.model}}}{Save models}
#'   \item{\code{\link{load.model}}}{Load models}}
#'
#' @docType package
#' @name NumberOne
NULL
#> NULL
