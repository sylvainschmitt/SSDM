#' SSDM: Stacked species distribution modelling
#'
#' SSDM is a package to map species richness and endemism based on stacked species distribution models (SSDM). Individual
#' SDMs can be created using a single or multiple algorithms (ensemble SDMs). For each species,
#' an SDM can yield a habitat suitability map, a binary map, a between-algorithm variance map,
#' and can assess variable importance, algorithm accuracy, and between-algorithm
#' correlation. Methods to stack individual SDMs include summing individual
#' probabilities and thresholding then summing. Thresholding can be based on a
#' specific evaluation metric or by drawing repeatedly from a Bernouilli
#' distribution. The SSDM package also provides a user-friendly interface
#' (\code{\link{gui}}).
#'
#' SSDM provides five categories of functions (that you can find in details
#' below): Data preparation, Modelling main functions, Model main methods, Model
#' classes, and Miscellaneous.
#'
#' @section Data preparation: \describe{\item{\code{\link{load_occ}}}{Load
#'   occurrence data} \item{\code{\link{load_var}}}{Load environmental variables}}
#'
#' @section Modelling main  functions: \describe{
#'   \item{\code{\link{modelling}}}{Build an SDM using a single algorithm}
#'   \item{\code{\link{ensemble_modelling}}}{Build an SDM that assembles
#'   multiple algorithms} \item{\code{\link{stack_modelling}}}{Build an SSDMs
#'   that assembles multiple algorithms and species}}
#'
#' @section Model main methods: \describe{
#'   \item{\code{\link{ensemble,Algorithm.SDM-method}}}{Build an ensemble SDM}
#'   \item{\code{\link{stacking,Ensemble.SDM-method}}}{Build an SSDM}
#'   \item{\code{\link{update,Stacked.SDM-method}}}{Update a previous SSDM with
#'   new occurrence data}}
#'
#'
#' @section Model classes: \describe{
#'   \item{\code{\linkS4class{Algorithm.SDM}}}{S4 class to represent SDMs}
#'   \item{\code{\linkS4class{Ensemble.SDM}}}{S4 class to represent ensemble
#'   SDMs} \item{\code{\linkS4class{Stacked.SDM}}}{S4 class to represent SSDMs}}
#'
#' @section Miscellaneous: \describe{ \item{\code{\link{gui}}}{User-friendly
#'   interface for SSDM package} \item{\code{\link{plot.model}}}{Plot SDMs}
#'   \item{\code{\link{save.model}}}{Save SDMs}
#'   \item{\code{\link{load.model}}}{Load SDMs}}
#'
#' @docType package
#' @name SSDM
NULL
#> NULL
