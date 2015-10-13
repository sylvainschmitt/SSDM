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
TreatVar <- function (Env, Reso.adapt = T, Norm = T) {
  cat('Variables treatment \n')
  
  # Resolution
  if(Reso.adapt) {
    cat('   resolution adaptation \n')
    reso = max(res(Env@layers[[1]]))
    # Coarser resolution measurement
    for (i in 1:length(Env@layers)){
      reso = max((res(Env@layers[[i]])),reso)
    }
    # Refine all stack resolution
    res(Env) = reso
  }
  
  # Normalizing variables
  if (Norm) {
    cat('   normalizing continuous variables \n\n')
    for (i in 1:length(Env@layers)) {
      #For not categorical variables
      if (!Env[[i]]@data@isfactor) { 
        Env[[i]] = Env[[i]]/Env[[i]]@data@max
      }
    }
  }
  
  return(Env)
}

#### Main #### ----
# Available algorithms :
# Working algorithms'GLM','GAM','MAXENT','ANN','CTA','GBM','RF', 'FDA', 'MARS
Env = LoadEnv()
Occurences = LoadOcc(Reduc = T)
Env = TreatVar(Env)
setwd("/home/amap/Documents/R/Sans Biomod 2/")

model1 = Modelling('GLM', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model2 = Modelling('GAM', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model3 = Modelling('MARS', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model4 = Modelling('CTA', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model5 = Modelling('GBM', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model62 = Modelling('RF', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model7 = Modelling('MAXENT', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model8 = Modelling('ANN', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')
model9 = Modelling('SVM', Occurences, Env, Xcol = 'Longitude', Ycol = 'Latitude')

obj = Algorithm.Niche.Model('GLM', data = data.frame(X = Occurences$Longitude, Y = Occurences$Latitude))
obj = PA.select(obj, Env)
obj = data.values(obj, Env, na.rm = T)
summary(obj@data)
obj = project(obj, Env)
obj = evaluate(obj)
obj = evaluate.axes(obj, thresh, Env)
