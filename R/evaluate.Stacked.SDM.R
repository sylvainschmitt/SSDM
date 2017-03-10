#' @include Stacked.SDM.R
#' @import methods
#' @importFrom raster reclassify rasterize extract stack
#' @importFrom sp SpatialPoints coordinates
NULL

#' Evaluate
#'
#' Evalaution of SSDM floristic composition with Pottier et al, 2013 method (see reference below)
#'
#' @param obj Stacked.SDM. SSDM to evaluate
#'
#' @return SSDM evaluation in a data.frame
#'
#' @name evaluate.Stacked.SDM
#'
#' @export
#'
#' @examples
#'
#' ## Not run:
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Evaluation
#' evaluate(SSDM)
#'
#' ## End(Not run)
#'
#' @references Pottier, J., Dubuis, A., Pellissier, L., Maiorano, L., Rossier, L., Randin, C.
#' F., … Guisan, A. (2013). The accuracy of plant assemblage prediction from
#' species distribution models varies along environmental gradients. Global
#' Ecology and Biogeography, 22(1), 52–63.
#' https://doi.org/10.1111/j.1466-8238.2012.00790.x
#'
#'
setMethod("evaluate", "Stacked.SDM", function(obj, ...){

  # Observed composition
  obs <- lapply(obj@enms, function(x){
    x <- x@data[x@data$Presence == 1,1:2]
  })
  obs <- mapply(function(x, n){
    x['species'] <- n ; return(x)
  }, n = lapply(strsplit(names(obj@enms), '.', fixed = TRUE), function(x){x[[1]]}),
  x= obs, SIMPLIFY = FALSE)
  obs <- do.call('rbind', obs)
  obs <- with(obs, table(paste(X, Y), species))
  obs <- as.data.frame.matrix(obs)
  obs[obs > 1] <- 1
  XY <- do.call('rbind', strsplit(row.names(obs), ' '))
  row.names(obs) <- 1:dim(obs)[1]
  XY <- apply(XY, 1, as.character)
  XY <- apply(XY, 1, as.numeric)
  XY <- data.frame(XY)
  names(XY) <- c('X','Y')

  # Predicted composition
  pred <- stack(lapply(obj@enms, function(x){
    x@binary
  }))
  names(pred) <- unlist(lapply(strsplit(names(obj@enms), '.', fixed = TRUE), function(x){x[[1]]}))
  pred <- extract(pred, XY)

  # Confusion matrix
  conftab <- obs*10 + pred
  conf <- data.frame(
    TP = apply(conftab, 1, function(x){sum(length(which(x == 11)))}), # True positive
    FN = apply(conftab, 1, function(x){sum(length(which(x == 10)))}), # False negative
    FP = apply(conftab, 1, function(x){sum(length(which(x == 1)))}), # False positive
    TN = apply(conftab, 1, function(x){sum(length(which(x == 0)))}) # True negative
  )

  # Evaluation indices
  n <- dim(conf)[2] # Renmaing to follow Pottier et al, 2013
  a <- conf$TP
  b <- conf$FP
  c <- conf$FN
  d <- conf$TN
  kappa.inter <- ((a + b)*(a + c) + (c + d)*(d + b)) / n*n
  eval <- data.frame(
    'species richness error' = rowSums(pred) - rowSums(obs),
    'prediction success' = (a + d) / n,
    'kappa' = (((a + d) / n) - kappa.inter) / (1 - kappa.inter),
    'specificity' = d / (b + d),
    'sensitivity' = a / (a + c),
    'Jaccard' = a / (a + b + c)
  )

  # Result
  eval <- data.frame(
    mean = apply(eval, 2, mean, na.rm = TRUE),
    SD = apply(eval, 2, sd, na.rm = TRUE)
  )

  return(data.frame(t(eval)))
})
