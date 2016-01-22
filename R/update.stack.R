#'@include ensemble_modelling.R stacking.R checkargs.R
NULL

#'Update a previous SSDM
#'
#'Update a previous SSDM with new occurrence data. The function takes as inputs
#'updated or new occurrence data from one species, previous environmental
#'variables, and an S4 \linkS4class{Stacked.SDM} class object containing a
#'previously built SSDM.
#'
#'@param object Stacked.SDM. The previously built SSDM.
#'@param Occurrences data frame. New or updated occurrence table (can be
#'  processed first by \code{\link{load_occ}}).
#'@param Env raster object. Environment raster object (can be processed first by
#'  \code{\link{load_var}}).
#'@param Xcol character. Name of the column  in the occurrence table  containing
#'  Latitude or X coordinates.
#'@param Ycol character. Name of the column in the occurrence table  containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence, by setting presence to 1 and
#'  absence to 0. If NULL presence-only dataset is assumed.
#'@param Spname character. Name of the new or updated species.
#'@param name character. Optional name given to the final SSDM produced, by
#'  default it's the name of the previous SSDM.
#'@param save logical. If set to true, the model is automatically saved.
#'@param path character. Name of the path to the directory to contain the saved
#'  SSDM.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'@param tmp logical. If set to true, the habitat suitability map of each
#'  algorithm is saved in a temporary file to release memory. But beware: if you
#'  close R, temporary files will be destroyed. To avoid any loss you can save
#'  your model with \code{\link{save.model}}.
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'@param ... additional parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'
#'@seealso \code{\link{stack_modelling}} to build SSDMs.
#'
#' @examples
#' \dontrun{
#' update(stack, Occurrences, Env, Spname = 'NewSpecie')
#' }
#'
#'@export
setMethod('update', 'Stacked.SDM',
          function(object,
                   # Modelling data input
                   Occurrences, Env,
                   # Occurrences reading
                   Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL, Spname = NULL,
                   # Model creation
                   name = stack@name, save = F, path = getwd(), thresh = 1001, tmp = F,
                   # Verbose
                   verbose = T, GUI = F,
                   # Modelling parameters
                   ...) {
            # Check arguments
            .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, Spname = Spname, name = name, save = save,
                         path = path, thresh = thresh, tmp = tmp, verbose = verbose, GUI = GUI)

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
                                     save = F, path = getwd(), PA, cv = stack@parameters$cv,
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
            enms['endemism'] = strsplit(stack@parameters$endemism, '|', fixed = 'T')[[1]]
            enms['rep.B'] = stack@parameters$rep.B
            newstack = do.call(stacking, enms)
            if(verbose) {cat('   done.\n')}

            if(!is.null(stack)) {
              # Paremeters
              newstack@parameters$sp.nb.origin = stack@parameters$sp.nb.origin

              # Saving
              if(save) {
                if(verbose){cat('Saving...\n')}
                if (!is.null(name)) {save.stack(newstack, name = name, path = path)}
                else {save.stack(newstack, path = path)}
                if(verbose) {cat('   done.\n')}
              }
            }

            return(newstack)
          })
