#' @include checkargs.R
#' @importFrom raster raster stack res extent crop reclassify as.factor
NULL

#'Load environmental variables data
#'
#'Function to load environmental variables data from rasters to do
#'\code{\link{Modelling}}, \code{\link{Ensemble.Modelling}} or
#'\code{\link{Stack.Modelling}}.
#'
#'@param directory character. Directory where are the environmental data
#'  variables files
#'@param files character. Files containing the environmental data variables, if
#'  NULL (default) all files with the precised format present in the directory
#'  will be loaded.
#'@param format character. Format used to load environmental data variables
#'  files (among .grd, .tif, .asc, .sdat, .rst, .nc, .tif, .envi, .bil, .img)
#'@param factors character. Environmental data variables which should be
#'  considered as factor variables
#'@param Norm logical. If true environmental data are normalized.
#'@param tmp logical. If true loaded environmental data variables rasters are
#'  read in temporary file avoiding to overload the random access memory. But
#'  beware, if you close R temporary files will be destructed and you'll need to
#'  reload your environmental datas directly from their files.
#'@param verbose logical. If true allow the function to print text in the
#'  console
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface) !
#'
#'@return A stack containing the environmental variables raster (normalized or
#'  not)
#'
#' @examples
#'\dontrun{
#' load.var(directory)
#'}
#'
#'@seealso \code{\link{load.occ}} to load occurences
#'
#'@export
load.var <- function (directory = getwd(), files = NULL,
                      format = c('.grd','.tif','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img'),
                      factors = NULL, Norm = T, tmp = T, verbose = T, GUI = F) {
  # Check arguments
  .checkargs(directory = directory, files = files, format = format, factors = factors,
             Norm = Norm, tmp = tmp, verbose = verbose, GUI = GUI)

  # pdir = getwd()
  if(verbose) {cat('Variables loading \n')}
  # setwd(directory)
  Env = stack()

  # Rasters loading
  files.null = files
  if (is.null(files)) {files.null = T} else {files.null = F}
  for (j in 1:length(format)) {
    if(files.null) {
      files = list.files(path = directory, pattern = paste0('.',format[j],'$'))
    }
    if (length(files) > 0) {
      for (i in 1:length(files)){
        Raster = raster(paste0(directory,'/',files[[i]]))
        # Extent and resolution check
        reso = res(Raster)
        extent = extent(Raster)
        if (j == 1  && i == 1) {
          resostack = reso
          extentstack = extent
        } else {
          resostack[1] = max(reso[1], resostack[1])
          resostack[2] = max(reso[2], resostack[2])
          # Extent and resolution adpatation
          extentstack@xmin = max(extentstack@xmin, extent@xmin)
          extentstack@xmax = min(extentstack@xmax, extent@xmax)
          extentstack@ymin = max(extentstack@ymin, extent@ymin)
          extentstack@ymax = min(extentstack@ymax, extent@ymax)
        }
        if(GUI) {incProgress(1/(length(files)*3), detail = paste(i,'loaded'))}
      }
    }
  }

  if(verbose) {cat('Variables treatment \n')}

  if(verbose) {cat('   resolution and extent adaptation...')}
  for (j in 1:length(format)) {
    if(files.null) {
      files = list.files(path = directory, pattern = paste0('.',format[j],'$'))
    }
    if (length(files) > 0) {
      for (i in 1:length(files)){
        Raster = raster(paste0(directory, '/', files[[i]]))
        Raster = reclassify(Raster, c(-Inf,-900,NA))
        Raster = crop(Raster, extentstack)
        names(Raster) = as.character(strsplit(files[i],paste0('.',format[j]))[1])
        if (names(Raster) %in% factors) {
          Raster = raster::as.factor(Raster)
          row.names(Raster@data@attributes[[1]]) = Raster@data@attributes[[1]]$ID
          fun = max
        } else {fun = mean}
        if (round(res(Raster)[1], digits = 6) != round(resostack[1], digits = 6) || round(res(Raster)[2], digits = 6) != round(resostack[2], digits = 6)) {
          cat(c((res(Raster)[1]/resostack[1]),(res(Raster)[2]/resostack[2])))
          Raster = aggregate(Raster, fact = (res(Raster)[1]/resostack[1]), fun = fun)
        }
        Env = stack(Env, Raster)
        if(GUI) {incProgress((1/(length(files)*3)), detail = paste(i,'treated'))}
      }
    }
  }
  if(verbose) {cat('   done... \n')}


  # Normalizing variables
  if (Norm) {
    if(verbose) {cat('   normalizing continuous variables \n\n')}
    for (i in 1:length(Env@layers)) {
      #For not categorical variables
      if (!Env[[i]]@data@isfactor) {
        Env[[i]] = Env[[i]]/Env[[i]]@data@max
      }
      if(GUI) {incProgress((1/length(Env@layers)/3), detail = paste(i,'normalized'))}
    }
  }

  # setwd(pdir)

  # Temporary files
  if (tmp) {
    path = get("tmpdir",envir=.PkgEnv)
    if (!("./.rasters" %in% list.dirs())) (dir.create(paste0(path,'/.rasters')))
    for (i in 1:length(Env@layers)) {Env[[i]] = writeRaster(Env[[i]], paste0(path,"/.rasters/", names(Env[[i]])), overwrite = T)}
  }

  return(Env)
}
