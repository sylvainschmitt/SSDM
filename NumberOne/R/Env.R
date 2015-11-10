#' Environmental variables ratser stack
#'
#' A raster stack containing 3 environmental variables on the scale of New
#' Caldeonia 'Grande Terre'. Climatic variables come from Bioclim datas, and
#' SUBSTRATE variable is from IRD atlas (see references below). To computing
#' time purpose rasters have been degraded to 5 kilometers resolution.
#'
#' @format A raster stack with 3 rasters: \describe{ \item{RAINFALL}{Annual mean
#'   rainfall} \item{TEMPERATURE}{Annual mean temperature} \item{SUBSTRATE}{Soil
#'   substrate (Categorical variable)} }
#'
#' @references Hijmans,  R.J.  &  Graham,  C.H.  (2006).  The  ability  of
#'   climate envelope  models  to  predict  the  effect  of  climate  change  on
#'   species distributions. Glob. Chang. Biol., 12, 1-10.
#'
#'   Fritsch, E. (2012) Les sols. Atlas de la Nouvelle-Caledonie (ed. by J.
#'   Bonvallot, J.-C. Gay and E. Habert), pp. 73-76. IRD-Congres de la
#'   Nouvelle-Caledonie, Marseille.
#'
"Env"
