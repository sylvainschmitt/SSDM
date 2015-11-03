#' @importFrom spThin thin
NULL

#' @export
load.occ = function(directory = getwd(), Env, file = NULL, ...,
                     Xcol = 'Longitude', Ycol = 'Latitude', Spcol = NULL,
                     GeoRes = T, reso = max(res(Env@layers[[1]])), verbose = T, GUI = F) {
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
