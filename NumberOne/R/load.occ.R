#' @include checkargs.R
#' @importFrom spThin thin
NULL

#'Load occurences data
#'
#'Function to load occurences data from a table to do \code{\link{Modelling}},
#'\code{\link{Ensemble.Modelling}} or \code{\link{Stack.Modelling}}.
#'
#'@param directory character. Directory where is the occurences table file
#'@param Env raster stack. Environment vraiables as a raster stack used to do
#'  the spatial thinning, it can be the result of the \code{\link{load.var}}
#'  function
#'@param file character. File containing the occurences table, if NULL (default)
#'  .csv file present in the directory will be loaded
#'@param ... additional parameters gave to \code{\link[utils]{read.csv}}
#'@param Xcol character. Name of the X coordinates column (Longitude)
#'@param Ycol character. Name of the Y coordinates column (Latitude)
#'@param Spcol character. Name of the column containing the species name or ID
#'@param GeoRes logical. Geographical thinning will be computed on occurences by
#'  species to avoid geographical bias in models
#'@param reso numeric. Resolution used to do the geographical thinnig, by
#'  default the environmental variables raster stack resolution (Env)
#'@param verbose logical. If true allow the function to print text in the
#'  console
#'@param GUI logical. Don't take that argument into account (parameter for the user interface) !
#'
#'@return A data frame containing the occurences data set (spatially thinned or
#'  not)
#'
#' @examples
#'\dontrun{
#' load.occ(directory)
#'}
#'
#'@seealso \code{\link{load.var}} to load environmental variables
#'
#'@export
load.occ = function(directory = getwd(), Env, file = NULL, ...,
                     Xcol = 'Longitude', Ycol = 'Latitude', Spcol = NULL,
                     GeoRes = T, reso = max(res(Env@layers[[1]])), verbose = T, GUI = F) {
  # Check arguments
  .checkargs(directory = directory, file = file, Xcol = Xcol, Ycol = Ycol, Spcol = Spcol,
             GeoRes = GeoRes, reso = reso, verbose = verbose, GUI = GUI)

  #pdir = getwd()
  if (verbose) {cat('Occurences loading \n')}
  #setwd(directory)
  if (is.null(file)) {
    file = as.character(list.files(path = directory, pattern = '.csv')[[1]])
  }
  file = paste0(directory,'/',file)
  Occurences = read.csv2(file = file, ...)  # Occ = occurences

  # Geographical resampling
  if (is.null(Spcol)) {
    Occurences$SpNULL = 1
    Spcol = 'SpNULL'
  }
  for (i in 1:length(levels(as.factor(Occurences[,which(names(Occurences)==Spcol)])))) {
    if (GeoRes) {
      if(verbose) {cat(levels(as.factor(Occurences[,which(names(Occurences)==Spcol)]))[i],'geographical resampling \n')}
      SpOccurences = subset(Occurences, Occurences[which(names(Occurences)==Spcol)] == levels(as.factor(Occurences[,which(names(Occurences)==Spcol)]))[i])
      thin.result = thin(SpOccurences, long.col = Xcol, lat.col = Ycol, spec.col = Spcol,
                         thin.par = reso, reps = 1, locs.thinned.list.return = T,
                         write.files = F, write.log.file = F, verbose = F)
      deleted = {}
      occ.indices = c(1:length(row.names(SpOccurences)))
      res.indices = as.numeric(row.names(thin.result[[1]]))
      for (i in 1:length(occ.indices)) {if(!(occ.indices[i] %in% res.indices)) {deleted = c(deleted, occ.indices[i])}}
      if (length(deleted) > 0) {Occurences = Occurences[-deleted,]}
    }
    if (Spcol == 'SpNULL') {Occurences = Occurences[-which(names(Occurences)=='SpNULL')]}
    if(GUI) {incProgress(1/length(levels(as.factor(Occurences[,which(names(Occurences)==Spcol)]))),
                            detail = paste(levels(as.factor(Occurences[,which(names(Occurences)==Spcol)]))[i],'thinned'))}
  }

  #setwd(pdir)
  return(Occurences)
}
