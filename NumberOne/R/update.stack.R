#'@include Ensemble.Modelling.R stacking.R checkargs.R
NULL

#'Update a previous stack
#'
#'Update a previous stack species ensemble distribution model with new
#'occurences datas. It takes in inputs an updated or new occurences data frame
#'from one specie, previous envrionmental variables, and an S4
#'\linkS4class{Stack.Species.Ensemble.Niche.Model} class object containing a
#'previously realized stack species model.
#'
#'@param object Stack.Species.Ensemble.Niche.Model. The previously realized stack
#'  species model
#'@param Occurences data frame. New or updated occurences table (can be treated
#'  first by \code{\link{load.occ}}).
#'@param Env raster object. Environnment raster object (can be treated first by
#'  \code{\link{load.var}}).
#'@param Xcol character. Name of the occurences table column containing Latitude
#'  or X coordinates.
#'@param Ycol character. Name of the occurences table column containing
#'  Longitude or Y coordinates.
#'@param Pcol character. Name of the occurences table column containing presence
#'  or absence value, if NULL presence-only data set is assumed.
#'@param Spname character. Name of the new or updated specie
#'@param name character. Optionnal name given to the final
#'  Stack.Specie.Ensemble.Niche.Model producted, by default it's the name of the
#'  previous stack.
#'@param save logical. If true the model is automatically saved.
#'@param directory character. If save is true, the name of the directory to save
#'  the model.
#'@param thresh numeric. binary map threshold computing precision parmeter, the
#'  higher it is the more accurate is the threshold but the longer is the
#'  modelling evaluation step !
#'@param tmp logical. If true ensemble models habitat suitability map and
#'  uncertainty map rasters are saved in temporary files to release memory.
#'@param verbose logical. If true allow the function to print text in the
#'  console
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface) !
#'@param ... additionnal parameters for the algorithm modelling function (see
#'  details below).
#'
#'@return an S4 \linkS4class{Stack.Species.Ensemble.Niche.Model} Class object
#'  viewable with \code{\link{plot.model}} method
#'
#'@seealso \code{\link{Stack.Modelling}} for stack species ensemble modelling
#'  with multiple algorithms and species
#'
#' @examples
#' \dontrun{
#' update(stack, Occurences, Env, Spname = 'NewSpecie')
#' }
#'
#'@export
setMethod(update, 'Stack.Species.Ensemble.Niche.Model',
          function(object,
                   # Modelling data input
                   Occurences, Env,
                   # Occurences reading
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
            ENM = Ensemble.Modelling(strsplit(stack@parameters$algorithms, '.', fixed = T)[[1]][-1],
                                     Occurences, Env, Xcol, Ycol, Pcol,
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
              i = which(names(stack@enms) == paste0(Spname,'.Ensemble.Niche.Model'))
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
