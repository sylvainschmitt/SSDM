#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(spThin)
library(biomod2)
library(rpart)
library(rgdal)
library(mgcv)
library(SDMTools)
library(e1071)
library(plotrix)
library(rgeos)
library(snowfall)

#### Data loading ####
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
    # Occurences = subset(Occurences,Occurences$Taxon == 'Zygogynum pancheri (Baill.) Vink')
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
    Env[[1]]@data@attributes[[1]]$ID = labels(Env[[1]]@data@attributes[[1]])[[1]]
    Env[[2]] = as.factor(Env[[2]])
  }
  return(Env)
}

#### Main ####
# Available algorithms :
# Working algorithms'GLM','GAM','MAXENT','ANN','CTA','GBM','RF', 'FDA', 'MARS
Env = LoadEnv()
Env = TreatVar(Env)
Occurences = LoadOcc(Reduc = T)
setwd("/home/amap/Documents/R/Essai5")

enm.list = sp.loop(Occurences, Env, models = 'GAM', PA = 'Min1', save = F, Norm = F, log = F)
proba.map(zoom = T)
proba.map(niche = T)
varimp()