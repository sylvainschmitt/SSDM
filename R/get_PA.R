#' @include Algorithm.SDM.R
#' @import methods
NULL

setGeneric("get_PA", function(obj) {
  return(standardGeneric("get_PA"))
})

setMethod("get_PA", "Algorithm.SDM", function(obj) {
  return(obj)
})

setMethod("get_PA", "GLM.SDM", function(obj) {
  PA <- list(nb = 1000, strat = "random")
  return(PA)
})

setMethod("get_PA", "GAM.SDM", function(obj) {
  PA <- list(nb = 1000, strat = "random")
  return(PA)
})

setMethod("get_PA", "MARS.SDM", function(obj) {
  PA <- list(nb = 1000, strat = "random")
  return(PA)
})

setMethod("get_PA", "CTA.SDM", function(obj) {
  PA <- list(nb = length(obj@data$Presence), strat = "random")
  return(PA)
})

setMethod("get_PA", "GBM.SDM", function(obj) {
  PA <- list(nb = length(obj@data$Presence), strat = "random")
  return(PA)
})

setMethod("get_PA", "RF.SDM", function(obj) {
  PA <- list(nb = length(obj@data$Presence), strat = "random")
  return(PA)
})

setMethod("get_PA", "MAXENT.SDM", function(obj) {
  PA <- list(nb = 10000, strat = "random")
  return(PA)
})

setMethod("get_PA", "ANN.SDM", function(obj) {
  PA <- list(nb = length(obj@data$Presence), strat = "random")
  return(PA)
})

setMethod("get_PA", "SVM.SDM", function(obj) {
  PA <- list(nb = length(obj@data$Presence), strat = "random")
  return(PA)
})
