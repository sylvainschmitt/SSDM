#### Results Toolbox ####

proba.map = function(list = enm.list, Occ = Occurences, nb = 1, niche = F, zoom = F) {
  i = 1
  title = 'Probability map'
  if (niche) {
    i = 2
    title = 'Niche map'
  }
  plot(list[[nb]]@proj[[i]], main = title)
  points(Occ$Longitude, Occ$Latitude, cex = 0.5, pch =16, col = 'red')
  if(zoom) {
    extent = drawExtent()
    plot(crop(list[[nb]]@proj[[i]], extent), main = title)
    points(Occ$Longitude, Occ$Latitude, cex = 0.5, pch =16, col = 'red')
  }
}

uncert = function(list = enm.list, Occ = Occurences, nb = 1, zoom = F) {
  plot(list[[nb]]@incert, main = 'Uncertaintiy map')
  points(Occ$Longitude, Occ$Latitude, cex = 0.5, pch =16, col = 'red')
  if(zoom) {
    extent = drawExtent()
    plot(crop(list[[nb]]@incert, extent), main = 'Uncertaintiy map')
    points(Occ$Longitude, Occ$Latitude, cex = 0.5, pch =16, col = 'red')
  }
}

varimp = function(list = enm.list, nb = 1) {
  barplot(list[[nb]]@axes.contrib$Mean, names.arg = row.names(list[[nb]]@axes.contrib), las = 2)
}

algo.perf = function(list = enm.list, nb = 1, metric = 'AUC') {
  # Available metrics : AUC, omission.rate, sensitivity, specificity, prop.correct, Kappa
  i = which(names(list[[nb]]@algo.eval) == metric)
  barplot(as.matrix(list[[nb]]@algo.eval[i]), 
          names.arg = names(list[[nb]]@algo.corr), 
          las = 2, main = metric)
}

load.enm = function(SpecieID) {
  directory = getwd()
  setwd(paste0(directory, "/",SpecieID,'_Results'))
  enm = ENM(
    proj = stack(raster('Rasters/Probability.tif'), raster('Rasters/Niche.tif')),
    incert = raster('Rasters/Uncertainity.tif'),
    eval = read.csv('Tables/ENMeval'),
    algo.eval = read.csv('Tables/AlgoEval'),
    algo.corr = read.csv('Tables/AlgoCorr'),
    axes.contrib = read.csv('Tables/VarImp')
  )
  setwd(directory)
  return(enm)
}

diversity = function(enms, choice = 'Niche') {
  if (choice == 'Probability') {
    proj = enms[[1]]@proj[[1]]
    for (i in 2:length(enms)) {proj = proj + enms[[i]]@proj[[1]]}
  }
  if (choice == 'Niche') {
    proj = enms[[1]]@proj[[2]]
    for (i in 2:length(enms)) {proj = proj + enms[[i]]@proj[[2]]}
  }
  if (choice == 'Uncertainity') {
    proj = enms[[1]]@incert
    for (i in 2:length(enms)) {proj = proj + enms[[i]]@incert}
  }
  plot(proj, main = choice)
  return(proj)
}