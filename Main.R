#### Sources #### ----
source('Data.prep.R')
source('Classes.R')
source('Modelling.R')

#### Data loading #### ----
OpenEnv <- function(EnvNames, format) {
  # Open environement variables
  Env = stack()
  for (i in 1:length(EnvNames)){
    Raster = raster(paste0(EnvNames[i],".",format))
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
LoadEnv <- function(directory = "/home/amap/Documents/Données/ZYGOGYNUM_DATASET/VARIABLES_ENVIRONNEMENTALES",
                    names = c("alizes","altitude","cti","ensoleillement","pente","routes","secteurs","substrat"),
                    format = 'tif') {
  setwd(directory)
  Env = OpenEnv(names, format)
  return(Env)
}
Env = LoadEnv()
Occurences = LoadOcc(Reduc = T)
Env = treat.var(Env)
Occurences = treat.occ(Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
setwd("/home/amap/Documents/R/Sans Biomod 2/")

#### Main #### ----
#ENM = Ensemble.Modelling(c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'),
#                        Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude', name = 'Specie', rep = 10)
MAXENT = Modelling('MAXENT', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude', name = 'Specie')
GLM = Modelling('GLM', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude', name = 'Specie')
ENM = ensemble(MAXENT, GLM)
plot(ENM)

