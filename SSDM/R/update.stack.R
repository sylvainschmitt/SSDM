#'@include ensemble_modelling.R stacking.R checkargs.R
NULL

#'Update a previous SSDM
#'
#'Update a previous SSDM with new occurrences datas. It takes in inputs updated
#'or new occurrences data frame from one specie, previous envrionmental data,
#'and an S4 \linkS4class{Stacked.SDM} class object containing a previously
#'realized stacked SDM.
#'
#'@param object Stacked.SDM. The previously realized SSDM.
#'@param Occurrences data frame. New or updated occurrences table (can be
#'  treated first by \code{\link{load_occ}}).
#'@param Env raster object. Environnment raster object (can be treated first by
#'  \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurences table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurences table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurences table specifying
#'  wether a line is a presence or an absence. If NULL presence-only data set is
#'  assumed.
#'@param Spname character. Name of the new or updated specie.
#'@param name character. Optionnal name given to the final SSDM produced, by
#'  default it's the name of the previous stack.
#'@param save logical. If true the model is automatically saved.
#'@param directory character. Name of the directory to contain the saved SSDM.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 & 1. The higher it is the more accurate
#'  is the threshold but the longer is the modelling evaluation step (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param tmp logical. If true the habitat suitability map of each algorithms is
#'  saved in a temporary file to release memory. But beware: if you close R,
#'  temporary files will be destroyed. To avoid any loss you can save your model
#'  with \code{\link{save.model}}.
#'@param verbose logical. If true allow the function to print text in the
#'  console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface) !
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Stacked.SDM} Class object viewable with
#'  \code{\link{plot.model}} method
#'
#'@seealso \code{\link{stack_modelling}} for SSDMs building.
#'
#' @examples
#' \dontrun{
#' update(stack, Occurrences, Env, Spname = 'NewSpecie')
#' }
#'
#'@export
setMethod(update, 'Stacked.SDM',
          function(object,
                   # Modelling data input
                   Occurrences, Env,
                   # Occurrences reading
                   Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL, Spname = NULL,
                   # Model creation
                   name = stack@name, save = F, directory = getwd(), thresh = 1001, tmp = F,
                   # Verbose
                   verbose = T, GUI = F,
                   # Modelling parameters
                   ...) {
            # Check arguments
            .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, Spname = Spname, name = name, save = save,
                         directory = directory, thresh = thresh, tmp = tmp, verbose = verbose, GUI = GUI)

            stack = object
            # New ENM creation
            if(verbose){cat('New specie ensemble distribution model creation...\n')}
            if(stack@parameters$PA == 'default') {
              PA = NULL
            } else {
              PA = list('nb'  = strsplit(stack@parameters$PA, '.', fixed = T)[[1]][1],
                        'strat' = strsplit(stack@parameters$PA, '.', fixed = T)[[1]][2])
            }
            if(!is.null(Spname)) {enm.name = Spname} else {enm.name = 'new_Specie'}
            ENM = ensemble_modelling(strsplit(stack@parameters$algorithms, '.', fixed = T)[[1]][-1],
                                     Occurrences, Env, Xcol, Ycol, Pcol,
                                     rep = as.numeric(stack@parameters$rep), enm.name,
                                     save = F, directory = getwd(), PA, cv = stack@parameters$cv,
                                     cv.param = as.numeric(strsplit(stack@parameters$cv.param, '|', fixed = T)[[1]][-1]),
                                     thresh = thresh, metric = stack@parameters$metric,
                                     axes.metric = stack@parameters$axes.metric,
                                     uncertainity = stack@uncertainity@data@haveminmax, tmp = tmp,
                                     ensemble.metric = strsplit(stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1],
                                     ensemble.thresh = as.numeric(strsplit(stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1]),
                                     weight = as.logical(stack@parameters$weight), ...)
            if(verbose) {cat('   done.\n')}

            # Test for new
            if(verbose){cat('Check if the specie already exist...\n')}
            if(!is.null(Spname)){
              i = which(names(stack@enms) == paste0(Spname,'.Ensemble.SDM'))
              if(verbose){cat(Spname,'replacement\n')}
              if(length(i) > 0) {
                stack@enms[[i]] = NULL
              } else {
                stack@parameters$sp.nb.origin = stack@parameters$sp.nb.origin + 1
              }
            } else {
              stack@parameters$sp.nb.origin = stack@parameters$sp.nb.origin + 1
            }
            if(verbose) {cat('   done.\n')}

            # New stacking
            if(verbose){cat('New stacking...\n')}
            enms = list()
            for(i in 1:length(stack@enms)) {enms[[i]] = stack@enms[[i]]}
            enms['method'] = stack@parameters$method
            enms['metric'] = stack@parameters$metric
            enms['thresh'] = thresh
            enms['rep.B'] = stack@parameters$rep.B
            newstack = do.call(stacking, enms)
            if(verbose) {cat('   done.\n')}

            if(!is.null(stack)) {
              # Paremeters
              newstack@parameters$sp.nb.origin = stack@parameters$sp.nb.origin

              # Saving
              if(save) {
                if(verbose){cat('Saving...\n')}
                if (!is.null(name)) {save.stack(newstack, name = name, directory = directory)}
                else {save.stack(newstack, directory = directory)}
                if(verbose) {cat('   done.\n')}
              }
            }

            return(newstack)
          })
