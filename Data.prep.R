##### Mask creation ##### ----
env.mask = function(Env) {
  mask = reclassify(Env[[1]], c(-Inf,Inf,0))
  for (i in 1:length(Env@layers)) {
    mask = mask + reclassify(Env[[i]], c(-Inf,Inf,1))
  }
  mask = mask / length(Env@layers)
  mask = reclassify(mask, c(-Inf,1,NA), right = F)
  return(mask)
}