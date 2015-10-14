#### Sources #### ----
source('Data.prep.R')
source('Classes.R')
source('Modelling.R')

#### Data loading #### ----
OpenEnv <- function(EnvNames) {
  # Open environement variables
  Env = stack()
  for (i in 1:length(EnvNames)){
    Raster = raster(paste0(EnvNames[i],".tif"))
    Raster = reclassify(Raster, c(-Inf,-900,NA))
    names(Raster) = EnvNames[i]
    
    Env = stack(Raster,Env)
  }
  return(Env)
}
LoadOcc <- function(Reduc = F) {
  setwd("/home/amap/Documents/Données/")
  Occurences = read.csv2(file = "ZYGOGYNUM_DATASET/OCCURENCES/occurences.csv", dec = ".")  # Occ = occurences
  Occurences$SpeciesID = as.factor(Occurences$SpeciesID)
  if(Reduc) {
    # Reducing to one specie for testing purposes
    Occurences = subset(Occurences,Occurences$SpeciesID == '6318')
    #Occurences = subset(Occurences,Occurences$Taxon == 'Zygogynum pancheri (Baill.) Vink')
    Occurences = droplevels(Occurences)
  }
  return(Occurences)
}
LoadEnv <- function(Reduc = F, Extent = extent(165.8179 , 167.3069 , -22.33978, -21.57512), Cat = T) {
  setwd("/home/amap/Documents/Données/ZYGOGYNUM_DATASET/VARIABLES_ENVIRONNEMENTALES")
  EnvNames = c("alizes","altitude","cti","ensoleillement","pente","routes","secteurs","substrat")
  if (Cat == F) {EnvNames = c("alizes","altitude","cti","ensoleillement","pente","routes")}
  Env = OpenEnv(EnvNames)
  if(Reduc) {
    # Decreasing resolution for testing purposes
    Env = crop(Env, extent(165.8179 , 167.3069 , -22.33978, -21.57512))
    Env = stack(Env)
  }
  if (Cat) {
    Env[[1]] = reclassify(Env[[1]], c(3,Inf,NA))
    names(Env[[1]]) = 'substrat'
    Env[[2]] = reclassify(Env[[2]], c(-Inf,0.5,NA))
    names(Env[[2]]) = 'secteur'
    Env[[1]] = as.factor(Env[[1]])
    Env[[2]] = as.factor(Env[[2]])
  }
  return(Env)
}
Env = LoadEnv()
Occurences = LoadOcc(Reduc = T)
Env = treat.var(Env)
Occurences = treat.occ(Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude', Spcol = 'SpeciesID')
setwd("/home/amap/Documents/R/Sans Biomod 2/")

#### Main #### ----
# Available algorithms : 'GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF', 'MAXENT', 'ANN', 'SVM'
# sink('Log', split = T)
# sink()
# print(enm)
# plot(enm)
# save.enm(enm)
# load.enm(enm)

enm = Ensemble.Modelling(c('GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF', 'MAXENT', 'ANN', 'SVM'),
                         Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude', rep = 1)

# Object handmaking ----
obj = Algorithm.Niche.Model('GBM', data = data.frame(X = Occurences$Longitude, Y = Occurences$Latitude))
obj@data$Presence = 1
obj@data$Train = F
obj@data$Train[sample.int(length(obj@data$Presence), round(length(obj@data$Presence)*0.7))] = T
obj = PA.select(obj, Env)
obj = data.values(obj, Env)
summary(obj@data)
obj = project(obj, Env)

