#' @include checkargs.R
#' @importFrom shiny incProgress
#' @importFrom spThin thin
#' @importFrom raster res extract
NULL

#'Load occurrence data
#'
#'Function to load occurrence data from a table to perform \code{\link{modelling}},
#'\code{\link{ensemble_modelling}} or \code{\link{stack_modelling}}.
#'
#'@param path character. Path to the directory that contains the occurrence table.
#'@param Env raster stack. Environmental variables in the form of a raster stack used to
#'  perform spatial thinning (can be the result of the
#'  \code{\link{load_var}} function).
#'@param file character. File containing the occurrence table, if NULL
#'  (default) the .csv file located in the path will be loaded.
#'@param ... additional parameters given to \code{\link[utils]{read.csv}}.
#'@param Xcol character. Name of the column  in the occurrence table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table  containing
#'  Longitude or Y coordinates.
#'@param Spcol character. Name of the column containing species names or IDs.
#'@param GeoRes logical. Geographical thinning will be perform on occurrences
#'  to limit geographical biases in the SDMs.
#'@param reso numeric. Resolution used to perform the geographical thinning, by
#'  default the resolution of the raster stack (Env).
#'@param verbose logical. If set to true, allows the function to print text in the
#'  console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'
#'@return A data frame containing the occurrence dataset (spatially thinned or
#'  not).
#'
#' @examples
#'\dontrun{
#' load.occ(path)
#'}
#'
#'@seealso \code{\link{load_var}} to load environmental variables.
#'
#'@export
load_occ = function(path = getwd(), Env, file = NULL, ...,
                     Xcol = 'Longitude', Ycol = 'Latitude', Spcol = NULL,
                     GeoRes = T, reso = max(res(Env@layers[[1]])), verbose = T, GUI = F) {
  # Check arguments
  .checkargs(path = path, file = file, Xcol = Xcol, Ycol = Ycol, Spcol = Spcol,
             GeoRes = GeoRes, reso = reso, verbose = verbose, GUI = GUI)

  #pdir = getwd()
  if (verbose) {cat('Occurrences loading \n')}
  #setwd(path)
  if (is.null(file)) {
    file = as.character(list.files(path = path, pattern = '.csv')[[1]])
  }
  if(!is.null(path)){file = paste0(path,'/',file)}
  Occurrences = read.csv2(file = file, ...)  # Occ = occurrences

  # Checking columns format
  if(!is.null(Spcol)){
    if(!inherits(Occurrences[,which(names(Occurrences) == Spcol)], 'factor')){
      Occurrences[,which(names(Occurrences) == Spcol)] = as.factor(Occurrences[,which(names(Occurrences) == Spcol)])
    }
  }
  if(!inherits(Occurrences[,which(names(Occurrences) == Xcol)], 'numeric')){
    if(inherits(Occurrences[,which(names(Occurrences) == Xcol)], 'factor')){
    Occurrences[,which(names(Occurrences) == Xcol)] = as.numeric(as.character(Occurrences[,which(names(Occurrences) == Xcol)]))
    Occurrences[,which(names(Occurrences) == Ycol)] = as.numeric(as.character(Occurrences[,which(names(Occurrences) == Ycol)]))
    }
  }

  # Checking points validity
  Occurrences$validity = raster::extract(Env[[1]], Occurrences[,c(which(names(Occurrences) == Xcol),which(names(Occurrences) == Ycol))])
  if(length(which(is.na(Occurrences$validity))) > 0){
    warning('You have occurrences that aren\'t in the extent of your environmental variables, they will be automatically removed ! \n')
    Occurrences = Occurrences[-which(is.na(Occurrences$validity)),]
  }
  Occurrences = Occurrences[-which(names(Occurrences) == 'validity')]
  Occurrences = droplevels(Occurrences)

  # Geographical resampling
  if (is.null(Spcol)) {
    Occurrences$SpNULL = 1
    Spcol = 'SpNULL'
    Occurrences$SpNULL = as.factor(Occurrences$SpNULL)
  }
  for (i in 1:length(levels(Occurrences[,which(names(Occurrences)==Spcol)]))) {
    if (GeoRes) {
      if(verbose) {cat(levels(as.factor(Occurrences[,which(names(Occurrences)==Spcol)]))[i],'geographical resampling \n')}
      SpOccurrences = subset(Occurrences, Occurrences[which(names(Occurrences)==Spcol)] == levels(as.factor(Occurrences[,which(names(Occurrences)==Spcol)]))[i])
      thin.result = thin(SpOccurrences, long.col = Xcol, lat.col = Ycol, spec.col = Spcol,
                         thin.par = reso, reps = 1, locs.thinned.list.return = T,
                         write.files = F, write.log.file = F, verbose = F)
      if(GUI) {incProgress(1/length(levels(as.factor(Occurrences[,which(names(Occurrences)==Spcol)]))),
                           detail = paste(levels(as.factor(Occurrences[,which(names(Occurrences)==Spcol)]))[i],'thinned'))}
      deleted = {}
      occ.indices = c(1:length(row.names(SpOccurrences)))
      res.indices = as.numeric(row.names(thin.result[[1]]))
      for (i in 1:length(occ.indices)) {if(!(occ.indices[i] %in% res.indices)) {deleted = c(deleted, occ.indices[i])}}
      deleted = row.names(SpOccurrences[deleted,])
      deleted = which(row.names(Occurrences) %in% deleted)
      if (length(deleted) > 0) {Occurrences = Occurrences[-deleted,]}
    }
  }
  Occurrences = droplevels(Occurrences)

  # Test species occurrences > 3
  for (i in 1:length(levels(Occurrences[,which(names(Occurrences)==Spcol)]))) {
    sp = levels(as.factor(Occurrences[,which(names(Occurrences)==Spcol)]))[i]
    spocc = subset(Occurrences, Occurrences[,which(names(Occurrences)==Spcol)] == sp)
    l = length(spocc[,1])
    if(l < 4){
      warning(paste(sp, 'have 3 or less occurrences after spatial thinning, it\'s not enough for modelling, this species will be automatically removed ! \n'))
      Occurrences = Occurrences[-which(Occurrences[,which(names(Occurrences)==Spcol)] == sp),]
      }
  }
  if (Spcol == 'SpNULL') {Occurrences = Occurrences[-which(names(Occurrences)=='SpNULL')]}
  Occurrences = droplevels(Occurrences)

  #setwd(pdir)
  return(Occurrences)
}
