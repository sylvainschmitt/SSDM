#' @include Ensemble.SDM.R checkargs.R
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPoints bbox
#' @importFrom raster raster stack reclassify mask calc overlay values rasterize rasterToPoints
NULL

#'Stack different ensemble SDMs in an SSDM
#'
#'This is a function to stack several ensemble SDMs in an SSDM. The function
#'takes as inputs several S4 \linkS4class{Ensemble.SDM} class objects produced
#'with \code{\link{ensemble_modelling}} or \code{\link{ensemble}} functions. The
#'function returns an S4 \linkS4class{Stacked.SDM} class object containing the
#'local species richness map, the between-algorithm variance map, and all
#'evaluation tables coming with (model evaluation, algorithm evaluation,
#'algorithm correlation matrix and variable importance), and a list of ensemble
#'SDMs for each species (see \code{\link{ensemble_modelling}}).
#'
#'@param enm,... character. Ensemble SDMs to be stacked.
#'@param name character. Optional name given to the final SSDM produced (by
#'  default 'Species.SDM').
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random bernoulli (\strong{B}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'@param richness df. If the method used to create the local species richness is
#'  among maximum likelyhood (\strong{ML}), probability ranking (\strong{PR}),
#'  trait range (\strong{TR}), or checkerboard (\strong{CB}) , richness
#'  parameter allow to input an observed richness data frame. Default is null
#'  implying a richness computed from occurence data.
#'@param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}). Needed only for stacking
#'  method using SESAM framework: probability ranking (\strong{PR}), trait range
#'  (\strong{TR}), or checkerboard (\strong{CB}).
#'@param range integer. Set a value of range restriction (in pixels) around
#'  presences occurrences on habitat suitability maps (all further points will
#'  have a null probability, see Crisp et al (2011) in references). If NULL, no
#'  range restriction will be applied.
#'@param endemism character. Define the method used to create an endemism map
#'  (see details below).
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'@param GUI logical. Don't take that argument into account (parameter for the
#'  user interface).
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'
#'  \strong{Methods:} Choice of the method used to compute the local species
#'  richness map (see Calabrez et al. (2014) for more informations, see
#'  reference below): \describe{\item{P}{(Probablity) sum probabilities of
#'  habitat suitability maps }\item{B}{(Random bernoulli) draw repeatedly from a
#'  Bernoulli distribution}\item{T}{(Threshold) sum the binary map obtained with
#'  the thresholding (depending on the metric, see metric
#'  parameter).}\item{ML}{(Maximum likelyhood) To
#'  describe}\item{PR}{(Probability ranking) To describe}\item{TR}{(Trait range)
#'  To describe}\item{CB}{(Checkerboard) To describe}}
#'
#'  \strong{Endemism:} Choice of the method used to compute the endemism map
#'  (see Crisp et al. (2001) for more information, see reference below):
#'  \describe{\item{NULL}{No endemism map}\item{WEI}{(Weighted Endemism Index)
#'  Endemism map built by counting all species in each cell and weighting each
#'  by the inverse of its range} \item{CWEI}{(Corrected Weighted Endemism Index)
#'  Endemism map built by dividing the weighted endemism index by the total
#'  count of species in the cell.}}First string of the character is the method
#'  either WEI or CWEI, and in those cases second string of the vector is used
#'  to precise range calculation, whether the total number of occurrences
#'  \strong{'NbOcc'} whether the surface of the binary map species distribution
#'  \strong{'Binary'}.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occ1 <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#' Occ2 <- subset(Occurrences, Occurrences$SPECIES == 'gracilis')
#'
#' # SSDM building
#' ESDM1 <- ensemble_modelling(c('CTA', 'SVM'), Occ1, Env, rep = 1,
#'                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                            name = 'elliptica', ensemble.thresh = c(0.6))
#' ESDM2 <- ensemble_modelling(c('CTA', 'SVM'), Occ2, Env, rep = 1,
#'                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                            name = 'gracilis', ensemble.thresh = c(0.6))
#' SSDM <- stacking(ESDM1, ESDM2)
#'
#' # Results plotting
#' plot(SSDM)
#' }
#'
#'@seealso \code{\link{stack_modelling}} to build SSDMs.
#'
#'@references C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005)
#'  "Selecting thresholds of occurrence in the prediction of species
#'  distributions." \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'  J.M. Calabrese, G. Certain, C.  Kraan, & C.F. Dormann (2014) "Stacking
#'  species distribution  models  and  adjusting  bias  by linking them to
#'  macroecological models." \emph{Global Ecology and Biogeography} 23:99-112
#'  \url{http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/calabrese2013globalecolbiogeogr.pdf}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'  M. D. Crisp, S. Laffan, H. P. Linder & A. Monro (2001) "Endemism in the
#'  Australian flora"  \emph{Journal of Biogeography} 28:183-198
#'  \url{http://biology-assets.anu.edu.au/hosted_sites/Crisp/pdfs/Crisp2001_endemism.pdf}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@rdname stacking
#'@export
setGeneric('stacking', function(enm, ..., name = NULL, method = 'P', rep.B = 1000, richness = NULL, Env = NULL, range = NULL, endemism = c('WEI','Binary'), verbose = TRUE, GUI = FALSE) {return(standardGeneric('stacking'))})

#' @rdname stacking
#' @export
setMethod('stacking', 'Ensemble.SDM', function(enm, ..., name = NULL, method = 'P', rep.B = 1000,
                                               richness = NULL, Env = NULL,
                                               range = NULL, endemism = c('WEI','Binary'),
                                               verbose = TRUE, GUI = FALSE) {

  # Check arguments
  .checkargs(enm = enm, name = name, method = method, rep.B = rep.B, richness = richness, range = range, endemism= endemism,
             verbose = verbose, GUI = GUI)

  enms = list(enm, ...)
  if (length(enms) < 2) {stop('You neeed more than one ensemble SDM to do stackings')}
  names = c()
  for (i in seq_len(length(enms))) {if(enms[[i]]@name %in% names) {stop('Ensemble models can\'t have the same name, you need to rename one of ',enms[[i]]@name)} else {names = c(names, enms[[i]]@name)}}
  if(verbose) {cat('Stack creation... \n')}
  stack = Stacked.SDM(diversity.map = reclassify(enm@projection[[1]], c(-Inf,Inf,0)),
                      endemism.map = reclassify(enm@projection[[1]], c(-Inf,Inf,0)),
                      uncertainty = reclassify(enm@uncertainty, c(-Inf,Inf,NA)),
                      parameters = enm@parameters)

  # Name
  if(verbose) {cat('   naming...')}
  if (is.null(name)) {name = 'Species'}
  stack@name = paste0(name,'.SSDM')
  if(verbose) {cat(' done. \n')}

  # Range restriction
  if(verbose) {cat('   range restriction...')}
  if(!is.null(range)) {
    for(j in seq_len(length(enms))) {
      nbocc = length(as.factor(enms[[j]]@data$Presence[enms[[j]]@data$Presence==1])) / sum(enms[[j]]@algorithm.evaluation$kept.model)
      occ = enms[[j]]@data[1:nbocc,]
      occ = occ[which(occ$Presence == 1),1:2]
      circles = list()
      for (i in seq_len(length(occ[,1]))) {
        x = occ$X[i]
        y = occ$Y[i]
        pts = seq(0, 2 * pi, length.out = 100)
        # xy = cbind(x + range/60 * sin(pts), y + range/60 * cos(pts))
        res = res(stack@endemism.map)[1]
        xy = cbind(x + range * res * sin(pts), y + range * res * cos(pts))
        circle = Polygon(xy)
        circles[i] = circle
      }
      sc= SpatialPolygons(list(Polygons(circles, 'Circles')))
      enms[[j]]@projection = mask(enms[[j]]@projection, sc, updatevalue = 0)
#       thresh = enms[[j]]@evaluation$threshold
#       enms[[j]]@projection = reclassify(enms[[j]]@projection, c(-Inf,thresh,0, thresh,Inf,1))
    }
  }
  if(verbose) {cat(' done. \n')}

  # Diversity map
  if(verbose) {cat('   diversity mapping...')}
  # Useless datacheck to prevent bugs to remove after debugging
  for (i in seq_len(length(enms))) {
    if(!inherits(enms[[i]]@projection, 'RasterLayer')){
      if(verbose) {cat('Error', enms[[i]]@name, 'is not a raster but a', class(enms[[i]]@projection)[1], '.\nIt will be removed for the stacking')}
      enms[[i]] = NULL
    }
  }
  if (method == 'P') { # Individual probabilities sum (Calabrese et al, 2014)
    if(verbose) {cat('\n Local species richness coomputed by summing individual probabilities. \n')}
    for (i in seq_len(length(enms))) {stack@diversity.map = stack@diversity.map + enms[[i]]@projection}
  }
  if (method == 'T') { # Threshold and sum (Calabrese et al, 2014)
    if(verbose) {cat('\n Local species richness coomputed by thresholding and then summing. \n')}
    for (i in seq_len(length(enms))) {
      stack@diversity.map = stack@diversity.map +
        reclassify(enms[[i]]@projection,
                   c(-Inf,enms[[1]]@evaluation$threshold,0, enms[[i]]@evaluation$threshold,Inf,1))}
  }
  if (method == 'B') { # Random Bernoulli distribution (Calabrese et al, 2014)
    if(verbose) {cat('\n Local species richness coomputed by drawing repeatedly from a Bernoulli distribution. \n')}
    proba = stack()
    for (i in seq_len(length(enms))) {proba = stack(proba, enms[[i]]@projection)}
    diversity.map = calc(proba, fun = function(...) {
      x = c(...)
      x[is.na(x)] = 0
      return(rbinom(lengths(x), rep.B, x))},
      forcefun = TRUE)
    stack@diversity.map = sum(diversity.map) / length(enms) / rep.B
  }

  ##NEW STACKING METHODS IN DEVELOPMENT##
  if(method %in% c('LH','PR','TR','CB')){ # Richness map needed

    if(is.null(richness)){
      Richness <- reclassify(enm@projection, c(-Inf,Inf,0))
      for (i in seq_len(length(enms))) {
        Richness <- rasterize(SpatialPoints(enms[[i]]@data[1:2]), Richness, field = enms[[i]]@data$Presence, update = TRUE)
        Richness <- reclassify(Richness, c(-Inf,0,NA))
      }
    } else {
      stop('Richness external input is not yet implemented !')
    }
    # Some case can't be fitted !
    if(all(values(Richness) %in% c(0, 1, NA))){
      stop('Observed Richness is always equal to 1, modelled richness can\'t be adjusted !')
    }

    if(method == 'LH'){ # Maximum likelyhood (Calabrese et al, 2014)
      # stop("Maximum likelyhood stacking method from Calabrese et al, 2014 is not yet implemented in this version of the pacakge.")
      for (i in seq_len(length(enms))) {stack@diversity.map = stack@diversity.map + enms[[i]]@projection}
      SSDM_Richness <- values(mask(stack@diversity.map, Richness))
      SSDM_Richness <- SSDM_Richness[-which(is.na(SSDM_Richness))]
      Richness <- values(Richness)
      Richness <- Richness[-which(is.na(Richness))]
      fit <- lm(Richness ~ SSDM_Richness)
      a <- fit$coefficients[1]
      b <- fit$coefficients[2]
      stack@diversity.map <- a + b * stack@diversity.map
    }

    if(method %in% c('PR', 'TR', 'CB')){ # MEM needed

      # Compute the MEM ith same ensemble than the ESDMs
      Richness <- data.frame(rasterToPoints(Richness, fun = function(x){x>0}))
      MEM <- ensemble_modelling(algorithms = unlist(strsplit(enm@parameters$algorithms, '.', fixed = TRUE))[-1],
                                Occurrences = Richness, Env = Env, Xcol = 'x', Ycol = 'y', Pcol = 'layer',
                                rep = enm@parameters$rep, name = 'MEM', cv = enm@parameters$cv,
                                cv.param = as.numeric(unlist(strsplit(enm@parameters$cv.param, '|', fixed = TRUE))[-1]),
                                metric = enm@parameters$metric, axes.metric = enm@parameters$axes.metric,
                                ensemble.metric = unlist(strsplit(enm@parameters$ensemble.metric, '.', fixed = TRUE))[-1],
                                ensemble.thresh = as.numeric(unlist(strsplit(enm@parameters$ensemble.thresh, '|', fixed = TRUE))[-1]),
                                uncertainty = FALSE, weight = as.logical(enm@parameters$weight), verbose = FALSE)

      if(method == 'PR'){ # Probability ranking (SESAM, D'Amen et al, 2015)
        # stop("Probability ranking stacking method from Calabrese et al, 2014 is not yet implemented in this version of the pacakge.")
        # Rank the species by probability
        # Sum species probability while Richness < RichnessMEM
        stack@diversity.map <- MEM@projection
      }

      if(method == 'TR'){ # Trait range (SESAM, D'Amen et al, 2015)
        stop("Trait range stacking method from D'Amen et al, 2015 is not yet implemented in this version of the pacakge.")
      }

      if(method == 'CB'){ # TChecerboard (SESAM, D'Amen et al, 2015)
        stop("Checkerboard stacking method from D'amen et al, 2015 is not yet implemented in this version of the pacakge.")
      }
    }
  }
  ##NEW STACKING METHODS IN DEVELOPMENT##

  names(stack@diversity.map) = 'diversity'
  if(verbose) {cat(' done. \n')}

  # uncertainty map
  if(verbose) {cat('   uncertainty mapping...')}
  uncertainities = stack()
  for (i in seq_len(length(enms))) {
    a = try(enms[[i]]@uncertainty)
    if (inherits(a, 'try-error')) {
      if(verbose) {cat('Ensemble model',enms[[i]]@name,'uncertinity map not computed')}
    } else {
        b = try(stack(uncertainities, a))
        if(inherits(b, 'try-error')) {
          if(verbose) {cat('Ensemble model', enms[[i]]@name, ':', b)}
          } else {
            uncertainities = b
          }
        }
  }
  a = try(calc(uncertainities, mean))
  if (inherits(a, 'try-error')) {if(verbose) {cat('No uncertainty map to do uncertainty mapping')}
  } else {
    stack@uncertainty = a
    names(stack@uncertainty) = 'uncertainty'
  }
  if(verbose) {cat(' done. \n')}

  # endemism map
  if(verbose) {cat('   endemism mapping...')}
  if(is.null(endemism)) {
    'unactivated'
  } else {
    for (i in seq_len(length(enms))) {
      if(endemism[2] == 'NbOcc'){
        endweight = length(as.factor(enms[[i]]@data$Presence[enms[[i]]@data$Presence==1])) / sum(enms[[i]]@algorithm.evaluation$kept.model)
      } else if(endemism[2] == 'Binary') {
        endweight = sum(values(reclassify(enms[[i]]@projection, c(-Inf,enms[[i]]@evaluation$threshold,0,enms[[i]]@evaluation$threshold,Inf,1))), na.rm = TRUE)
      }
      if(endemism[1] == 'WEI') {
        stack@endemism.map = stack@endemism.map + enms[[i]]@projection / endweight
      } else if (endemism[1] == 'CWEI') {
        stack@endemism.map = stack@endemism.map + overlay(enms[[i]]@projection, stack@diversity.map, fun =
                                                            function(x,y){
                                                              y = round(y)
                                                              x[which(y > 0)] =  x[which(y > 0)] / endweight / y[which(y > 0)]
                                                              return(x)})
      }
    }
    stack@endemism.map = stack@endemism.map / stack@endemism.map@data@max
  }
  if(verbose) {cat(' done. \n')}

  # Evaluation
  if(verbose) {cat('   evaluating...')}
  stack@evaluation = enm@evaluation
  for (i in 2:length(enms)) {stack@evaluation = rbind(stack@evaluation, enms[[i]]@evaluation)}
  a = stack@evaluation[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in seq_len(length(stack@evaluation))) {a[i] = c(mean(stack@evaluation[,i], na.rm = TRUE), sd(stack@evaluation[,i], na.rm = TRUE))}
  stack@evaluation = a
  if(verbose) {cat(' done. \n')}

  # variable Importance
  if(verbose) {cat('   comparing variable importnace...')}
  stack@variable.importance = enm@variable.importance
  for (i in 2:length(enms)) {
    a  = try(rbind(stack@variable.importance, enms[[i]]@variable.importance))
    if (inherits(a, 'try-error')) {cat(a)} else {stack@variable.importance = a}
  }
  a = stack@variable.importance[1:2,]
  row.names(a) = c('Mean', 'SD')
  for (i in seq_len(length(stack@variable.importance))) {a[i] = c(mean(stack@variable.importance[,i]), sd(stack@variable.importance[,i]))}
  stack@variable.importance = a
  if(verbose) {cat(' done. \n')}

  # Algorithm Correlation
  if(verbose) {cat('   comparing algorithms correlation...')}
  algo = c() # Listing all algorithms presents in enms and renaming enms row and columns
  for (i in seq_len(length(enms))) {
    if(length(enms[[i]]@algorithm.correlation) == 0) {if(verbose) {cat('\n', enms[[i]]@name,'algorithms correlation has not been computed. \n')}} else {
      for (j in seq_len(length(enms[[i]]@algorithm.correlation))) {
        if (length(strsplit(names(enms[[i]]@algorithm.correlation)[j], '.', fixed = TRUE)[[1]]) > 1){
          a = strsplit(row.names(enms[[i]]@algorithm.correlation)[j], '.SDM', fixed = TRUE)[[1]][1]
          a = tail(strsplit(a, '.', fixed = TRUE)[[1]], n = 1)
          names(enms[[i]]@algorithm.correlation)[j] = a
          row.names(enms[[i]]@algorithm.correlation)[j] = a
        }
        if(!(names(enms[[i]]@algorithm.correlation)[j] %in% algo)) {algo = c(algo, names(enms[[i]]@algorithm.correlation)[j])}
      }
    }
  }
  mcorr = data.frame(matrix(nrow = length(algo), ncol = length(algo)))
  names(mcorr) = algo
  row.names(mcorr) = algo
  if(length(algo) > 0) {
    for (i in seq_len(length(algo))) {
      for (j in seq_len(length(algo))) {
        if(i > j) {
          corr = c()
          for (k in seq_len(length(enms))) {
            if(length(enms[[k]]@algorithm.correlation) != 0) {
              row = which(row.names(enms[[k]]@algorithm.correlation) == row.names(mcorr)[j])
              col = which(names(enms[[k]]@algorithm.correlation) == names(mcorr)[i])
              if(length(row) > 0 && length(col) > 0) {corr = c(corr, enms[[k]]@algorithm.correlation[row,col])}
            }
            mcorr[i,j] = mean(corr, na.rm = TRUE)
          }
        }
      }
    }
  }
  stack@algorithm.correlation = mcorr
  if(verbose) {cat(' done. \n')}

  # Algorithm Evaluation
  if(verbose) {cat('   comparing algorithms evaluation')}
  stack@algorithm.evaluation = enm@algorithm.evaluation
  for (i in 2:length(enms)) {stack@algorithm.evaluation = rbind(stack@algorithm.evaluation, enms[[i]]@algorithm.evaluation)}
  stack@algorithm.evaluation$algo = 'algo'
  for (i in seq_len(length(row.names(stack@algorithm.evaluation)))) {stack@algorithm.evaluation$algo[i] = strsplit(row.names(stack@algorithm.evaluation),'.', fixed = TRUE)[[i]][2]}
  stack@algorithm.evaluation = aggregate.data.frame(
    stack@algorithm.evaluation[-9],
    by = list(stack@algorithm.evaluation[,which(names(stack@algorithm.evaluation) == 'algo')]),
    FUN = mean)
  row.names(stack@algorithm.evaluation) = stack@algorithm.evaluation$Group.1
  stack@algorithm.evaluation = stack@algorithm.evaluation[-1]
  stack@algorithm.evaluation = stack@algorithm.evaluation[-which(names(stack@algorithm.evaluation) == 'algo')]
  if(verbose) {cat(' done. \n')}

  # ENMS
  for (i in seq_len(length(enms))) {suppressWarnings({stack@enms[enms[[i]]@name] = enms[[i]]})}

  # Parameters
  stack@parameters$method = method
  if (method == 'B') {stack@parameters$rep.B = rep.B}
  stack@parameters$range = range
  stack@parameters$endemism = paste0(endemism[1],'|',endemism[2])

  return(stack)})
