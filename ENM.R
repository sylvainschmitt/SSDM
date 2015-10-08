#### ENM class ####
setClass('ENM', 
         representation(name = 'character', proj = 'RasterStack', eval = 'data.frame', data = 'data.frame', incert = 'RasterLayer', algo.eval = 'data.frame', algo.corr = 'data.frame', axes.contrib = 'data.frame'), 
         prototype(name = character(), proj = stack(raster()), eval= data.frame(), data = data.frame(), incert = raster(), algo.eval = data.frame(), algo.corr = data.frame(), axes.contrib = data.frame()))
ENM <- function(name = character(), proj = stack(raster()), incert = raster(), eval = data.frame(), algo.eval = data.frame(), algo.corr = data.frame(), axes.contrib = data.frame(), data = data.frame()) {
  return(new('ENM', name = name, proj = proj, incert = incert, eval = eval, data = data,  algo.eval = algo.eval, algo.corr = algo.corr, axes.contrib = axes.contrib))
}

# 1 - Occurences treatment #
TreatOcc <- function(Occurences, 
                     Spcol = 'SpeciesID', 
                     Xcol = 'Longitude', 
                     Ycol = 'Latitude'
                     ) {
  cat('Occurences treatment \n\n')
  i = which(names(Occurences) == Spcol)
  j = which(names(Occurences) == Xcol)
  k = which(names(Occurences) == Ycol)
  Occurences = as.data.frame(cbind(Occurences[c(i, j, k)]))
  names(Occurences)[1] = 'Sp'
  names(Occurences)[2] = 'X'
  names(Occurences)[3] = 'Y'
  Occurences$Sp = as.factor(Occurences$Sp)
  return(Occurences)
}

# 2 -  Variables treatment #
# Without multiscaling for the moment
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

# 3 - Resampling #
# Geographical resampling with thin {spThin}
Resample = function(Occurences, Env, GeoRes = T, reso = max(res(Env@layers[[1]]))) {
  if (GeoRes) {
    cat('Geographical resampling \n\n')
    Occ_resamp = thin(Occurences, lat.col = "X", long.col = "Y", spec.col = "Sp", thin.par = reso, reps = 1, locs.thinned.list.return = T, write.files = F, write.log.file = F, verbose = F)[[1]]
    Occ_resamp$Sp = Occurences$Sp[1]
  }
  return(Occ_resamp)
}

# 4 - Modelling #

# Pseudo-absences and data preparation #
points.in.categories = function(Env, nb = 1000) {
  points = data.frame(matrix(nrow = 0, ncol = 2))
  names(points) = c('X', 'Y')
  for(c in 1:length(Env@layers)) {
    if(Env[[c]]@data@isfactor) {
      for(i in Env[[c]]@data@min:Env[[c]]@data@max) {
        point = data.frame(matrix(nrow = 0, ncol = 2))
        names(point) = c('X', 'Y')
        Secteur = reclassify(Env[[c]], c(-Inf,(i-0.5),NA, (i+0.5),Inf,NA))
        #border = rasterToPolygons(Secteur, dissolve = T)
        border = extent(Env[[1]])
        while(length(point$x)==0) {
          x = runif(nb, min = bbox(border)[1,1], max = bbox(border)[1,2])
          y = runif(nb, min = bbox(border)[2,1],max = bbox(border)[2,2])
          point = data.frame(x = x, y = y)
          coordinates(point) <- ~ x + y
          value = extract(Secteur, point)
          point = as.data.frame(point)
          point = point[-which(is.na(value)),]
        }
        point = point[1,]
        points = rbind(points,point)
      }
    }
  }
  return(points)
}

data.prep = function(PA.rep, PA.nb, Occurences, Env, PA.strat = 'random', PA.dist = NULL) {
  cat('Data preparation \n')
  
  # Mask defining
  cat('   mask defining')
  border = readOGR(dsn = 'bord', layer = 'bord', verbose = F)
  border = crop(border, extent(Env))
  if (PA.strat == 'random') {
    cat('   random selection \n')
  }
  if (PA.strat == '2nd') {
    cat('   second far selection \n')
    circles = list()
    for (i in 1:length(Occurences$Latitude)) {
      x = Occurences$Latitude[i]
      y = Occurences$Longitude[i]
      n = 100
      r = PA.dist/60
      pts = seq(0, 2 * pi, length.out = n)
      xy = cbind(x + r * sin(pts), y + r * cos(pts))
      circle = Polygon(xy)
      circles[i] = circle
    }
    sl= SpatialPolygons(list(Polygons(circles, 'Circles')))
    border = gIntersection(sl, border)
  }
  
  # Pseudo-Absences selection
  cat('   pseudo-absences selection : \n')
  data = as.data.frame(cbind(Occurences$Latitude, Occurences$Longitude))
  names(data) = c('x','y')
  PAcat = points.in.categories(Env) # Selection of one PA by each category
  data = rbind(data, PAcat)
  l = length(data$x)
  for (i in 1:PA.rep) {
    PA = data.frame(matrix(nrow = 0, ncol = 2))
    names(PA) = c('x','y')
    while (length(PA[,1]) < PA.nb) {
      x = runif(PA.nb*10, min = bbox(border)[1,1], max = bbox(border)[1,2])
      y = runif(PA.nb*10, min = bbox(border)[2,1],max = bbox(border)[2,2])
      points = data.frame(X = x, Y = y)
      coordinates(points) <- ~ X + Y
      check = over(points,border)
      if (PA.strat == 'random') {
        x = x[-which(is.na(check$GRIDCODE))]
        y = y[-which(is.na(check$GRIDCODE))]
      }
      if (PA.strat == '2nd') {
        x = x[-which(is.na(check))]
        y = y[-which(is.na(check))]
      }
      PA = rbind(PA, data.frame(x,y))
      cat((length(PA[,1])+PA.nb*i),' ')
    }
    PA = PA[1:PA.nb,]
    data = rbind(data, PA)
  }
  
  # Table creation
  cat('   table creation \n')
  data$Presence = 0
  data$Presence[1:length(Occurences$Longitude)] = 1
  for (i in 1:PA.rep) {
    data$Run = F
    data$Run[1:l] = T
    init = length(Occurences$Longitude) + (i-1)*PA.nb
    data$Run[init:(init+PA.nb)] = T
    names(data)[length(data)] = paste0('Run',i)
  }
  
  # Values extraction
  cat('   values extraction \n')
  data = cbind(data, extract(Env, cbind(data$x, data$y)))
  
  # Categorical variables as factor
  for (i in 1:length(Env@layers)) {
    if(Env[[i]]@data@isfactor) {
      col = which(names(data) == Env[[i]]@data@names)
      data[,col] = as.factor(data[,col])
    }
  }
  
  cat('done \n\n')
  return(data)
}


# Modeling #
Modelling = function (Occurences, Env, 
                      models = c('all'),
                      PA = 'Adapt', # Adapt, Min1, Min2, 10M, 10, Disk1
                      test = 'AIC',
                      epsilon = 1e-08,
                      maxit = 500,
                      cv = 3,
                      trees = 2500,
                      thresh.shrink = 1e-03,
                      final.leave = 1,
                      split = 70
                      ) {
  
  # Parameterization
  if(models == 'all') {models = c('GLM','GAM','MAXENT','ANN','CTA','GBM','RF','FDA','MARS')} # 'FDA' et 'MARS' en débugage, 'SVM' en rajout
  
  # Initalization
  models.out = list() #Results list
  sp = as.character(Occurences$Sp[1]) # Defining specie
  cat('Modeling step \n')
  
  # Model creation
  failed = c()
  for (i in 1:length(models)) {
    cat('   model',models[i],'...\n')
    if(PA == 'Adapt') {
      if(models[i] %in% c('GLM','GAM','MAXENT','ANN','SVM')){
        PA.rep = 10
        PA.nb = 100
        PA.strat = 'random'
        PA.dist = NULL
      }
      if(models[i] == 'MARS'){
        PA.rep = 10
        PA.nb = 100
        PA.strat = 'random'
        PA.dist = NULL
      }
      if(models[i] == 'FDA'){
        PA.rep = 10
        PA.nb = 100
        PA.strat = '2nd'
        PA.dist = 2
      }
      if(models[i] %in% c('CTA','GBM','RF')){
        if(length(Occurences[,1]) <10000) {
          PA.rep = 10
          PA.nb = length(Occurences[,1])
          PA.strat = '2nd'
          PA.dist = 2
        } else {
          PA.rep = 1
          PA.nb = length(Occurences[,1])
          PA.strat = '2nd'
          PA.dist = 2
        }
      }
    }
    if(PA == 'Min2') {
      PA.rep = 2
      PA.nb = 100
      PA.strat = 'random'
      PA.dist = NULL
    }
    if(PA == 'Min1') {
      PA.rep = 1
      PA.nb = 100
      PA.strat = 'random'
      PA.dist = NULL
    }
    if(PA == '10') {
      PA.rep = 1
      PA.nb = 10
      PA.strat = 'random'
      PA.dist = NULL
    }
    if(PA == '10M') {
      PA.rep = 1
      PA.nb = 10000
      PA.strat = 'random'
      PA.dist = NULL
    }
    if(PA == 'Disk1') {
      PA.rep = 1
      PA.nb = 100
      PA.strat = '2nd'
      PA.dist = 2
    }
    
    # Data preparation
    cat('       formating data \n')
    data = data.prep(PA.rep, PA.nb, Occurences = Occurences, Env = Env, PA.strat = PA.strat, PA.dist = PA.dist)
    # Check PAs and sectors
    plot(Env[[1]], main = 'Occurences and Pseudo-absences')
    points(data$x, data$y, pch = 16, col = as.factor(PA$Presence), cex = 0.5)
    
    # Modeling
    cat('       modeling \n')
    model = modeling()
    if(!inherits(model,"try-error")) {
      models.out[[i]] = list(sp, models[i], data, model)
    } else {
      cat('   ',models[i],'modeling failed \n')
      cat(model, '\n')
      failed = c(failed, i)
    }
    cat('   done \n')
    # models.out[[i]] [[1]] : Sp, [[2]] : Algo, [[3]] : Data, [[4]] : Model
  }
  cat('Modeling step done \n\n')
  
  # Removing failed models
  if (length(failed) > 0) {
    for (i in 1:length(failed)) {
      models.out[failed[i]] = NULL
      failed = failed - 1
    }
  }
  
 return(models.out)
}

# 5 - Projection creation #
Project = function (Models) {
  cat('Projecting ...\n')
  failed = c()
  for (i in 1:length(Models)) {
    cat(Models[[i]][[2]],'\n')
    Proj = try(BIOMOD_Projection(Models[[i]][[4]], Env, paste(Models[[i]][[1]], "Proj"), selected.models = 'all', binary.meth = 'ROC', build.clamping.mask = F, silent = T))
    if(!inherits(Proj,"try-error")) {
      Models[[i]][[5]] = Proj
    } else {
      cat('   ',Models[[i]][[2]],'projection failed \n')
      cat(Proj, '\n')
      failed = c(failed, i)
    }
    # Models[[i]] [[1]] : Sp, [[2]] : Algo, [[3]] : Data, [[4]] : Model, [[5]] : Projection
  }
  
  # Removing failed projections
  if (length(failed) > 0) {
    for (i in 1:length(failed)) {
      Models[failed[i]] = NULL
      failed = failed - 1
    }
  }
  
  cat('Projections done \n\n')
  return(Models)
}

# 6 - Algorithms models ensemble forecasting #
AlgoENM <- function(Models, AUC.thresh = 0.75, thresh = 1001) {
  # What the functions do :
  # Models[[i]] [[1]] : Sp, [[2]] : Algo, [[3]] : Data, [[4]] : Model, [[5]] : Projections * Algo repetitions
  # Transformed to
  # Models[[i]] [[1]] : Sp, [[2]] : Algo, [[3]] : Data, [[4]] : Model, [[5]] : Algorithm ensemble forecasting projection
  
  erased = c() # Vector for algorithms to erase
  
  for (i in 1:length(Models)) {
    Algo.ENM = ENM()
    
    cat(Models[[i]][[2]],'for specie',Models[[i]][[1]],'ensemble forecasting \n')
    
    # Evaluation
    Eval = get_evaluations(Models[[i]][[4]], as.data.frame = T)
    
    cat('  Projection forecasting \n')
    
    # Projections fusion
    Proj = Models[[i]][[5]]@proj@val
    ENM = reclassify(subset(Proj,1), c(-Inf,Inf,0)) # Null projection
    MaxVal = 0 # Null Proba value
    for (j in 1:length(Proj@layers)) {
      Nom_proj = Proj@layers[[j]]@data@names
      Nom_frac = strsplit(Nom_proj, split = "_")
      Nom_eval = paste0(Nom_frac[[1]][4],"_",Nom_frac[[1]][3],"_",Nom_frac[[1]][[2]])
      AUC = Eval[which(Eval$Model.name == Nom_eval),3]
      cat('    layers',j,'AUC =',AUC,':')
      if (AUC > AUC.thresh) {
        ENM = ENM + Proj@layers[[j]]*AUC
        MaxVal = MaxVal + AUC*1000
        cat(' accepted\n')
      } else {
        cat(' rejected\n')
      }
    }
    ENM = ENM/MaxVal # Scaling projection
    Algo.ENM@proj = stack(ENM)
    
    
    # Observation and Predictions
    Obs = Models[[i]][[3]]@data.species
    Obs[is.na(Obs)] = 0
    Coord = Models[[i]][[3]]@coord
    Pred = extract.data(Coord, ENM)
    Pred[is.na(Pred)] = 0
    Algo.ENM@data = cbind(Obs,Coord)
    
    # Threshold computing
    cat('  Threshold computing \n')
    seuil = optim.thresh(obs = Obs, pred = Pred, threshold = thresh)
    seuil = mean(seuil$`max.sensitivity+specificity`)
    
    # Projection evaluation
    cat('  Evaluation \n')
    Algo.ENM@eval = accuracy(Obs, Pred, seuil)
    
    # Variables importance
    cat('  Variables importance \n')
    VarImp = as.data.frame(Models[[i]][[4]]@variables.importances@val)
    Col = length(VarImp[1,])
    for (k in 1:length(VarImp[,1])) {
      x = {}
      for (l in 1:Col) {
        x = c(x, VarImp[k,l])
      }
      VarImp$Mean[k] = mean(x, na.rm = T)
    }
    Algo.ENM@axes.contrib = as.data.frame(VarImp[Col+1])
    
    # Saving or deleting the object if no layers kept
    if (Algo.ENM@proj[[1]]@data@min == Inf) {
      cat('  Algorithm erased due to not sufficient AUC \n \n')
      erased = c(erased,i)
    } else {
      Models[[i]][[5]] = Algo.ENM
      cat('  Done and saved \n \n')
    }
  }
  
  # Erasing algorithms rejected
  if (length(erased) > 0) {
    for (i in 1:length(erased)) {
      Models[[erased[i]]] = NULL
      erased = erased - 1
    }
  }
  
  return(Models)
}

# 7 - Ensemble forecasting #
ENMproj <- function(Models, thresh = 1001) {
  # Models[[i]] [[1]] : Sp, [[2]] : Algo, [[3]] : Data, [[4]] : Model, [[5]] : Projection (object of class ENM)
  
  enm = ENM()
  
  cat('Ensemble forecasting for algorithms : \n')
  
  # Projections, data, variable importances, and algorithm evaluation forecasting
  Proj = reclassify(Models[[1]][[5]]@proj[[1]], c(-Inf,Inf,0))
  MaxVal = 0
  data = data.frame(matrix(ncol = 3, nrow = 0))
  algo.eval = data.frame(matrix(ncol = 7, nrow = 0))
  VarImp = data.frame(matrix(ncol = 0, nrow = length(Models[[1]][[5]]@axes.contrib[,1])))
  row.names(VarImp) = row.names(Models[[1]][[5]]@axes.contrib)
  names(data) = names(Models[[1]][[5]]@data)
  names(algo.eval) = names(Models[[1]][[5]]@eval)
  for (i in 1:length(Models)) {
    cat('  ', Models[[i]][[2]], '...')
    AlgoENM = Models[[i]][[5]]
    Proj = Proj + AlgoENM@proj*AlgoENM@eval$AUC
    MaxVal = MaxVal + AlgoENM@eval$AUC
    data = rbind(data, AlgoENM@data)
    row.names(Models[[i]][[5]]@eval) = Models[[i]][[2]]
    algo.eval = rbind(algo.eval, AlgoENM@eval)
    names(Models[[i]][[5]]@axes.contrib) = Models[[i]][[2]]
    VarImp =cbind(VarImp, Models[[i]][[5]]@axes.contrib)
    cat(' done \n')
  }
  Proj = Proj / MaxVal
  enm@proj = stack(Proj)
  enm@data = data
  enm@algo.eval = algo.eval
  
  # Algorithms correlation
  cat('Algorithms correlation measurment  \n')
  Proj = stack()
  for (i in 1:length(Models)) {
    names(Models[[i]][[5]]@proj[[1]]) = Models[[i]][[2]]
    Proj = stack(Proj, Models[[i]][[5]]@proj[[1]])
  }
  enm@algo.corr = as.data.frame(layerStats(Proj, 'pearson', na.rm = T)$`pearson correlation coefficient`)
  
  
  # Observations and predictions
  Obs = data$Obs
  Obs[is.na(Obs)] = 0
  Coord = data[2:3]
  Pred = extract.data(Coord, enm@proj[[1]])
  Pred[is.na(Pred)] = 0
  
  # Threshold computing
  cat('Threshold computing \n')
  seuil = optim.thresh(obs = Obs, pred = Pred, threshold = thresh)
  seuil = mean(seuil$`max.sensitivity+specificity`)
  
  # Model evaluation
  cat('Model evaluation \n')
  enm@eval = accuracy(Obs, Pred, seuil)
  
  # Uncertainties map computing
  cat('Uncertainties computing \n')
  if (length(Proj@layers) > 1) {
    enm@incert = calc(Proj, var)
    names(enm@incert) = 'Uncertainity map'
  }
  
  # Presence/absence map computing
  cat('Binary transformation \n')
  enm@proj = stack(enm@proj, BinaryTransformation(enm@proj[[1]], enm@eval$threshold))
  names(enm@proj) = c('Probability map','Niche map')
  
  # Variables importances
  cat('Variables importances \n\n')
  Col = length(VarImp[1,])
  for (k in 1:length(VarImp[,1])) {
    x = {}
    for (l in 1:Col) {
      x = c(x, VarImp[k,l])
    }
    
    VarImp$Mean[k] = mean(x, na.rm = T)
    VarImp$Var[k] = var(x, na.rm = T)
  }
  enm@axes.contrib = as.data.frame(cbind(VarImp[Col+1],VarImp[Col+2]))
  
  return(enm)
}

# 8 - Saving #
save.enm = function(enm, directory = as.character(getwd()),
                    # Objects to save
                    Probabilities = T,
                    Niche = T,
                    Uncertainity = T,
                    ENMeval = T,
                    AlgoEval = T,
                    AlgoCorr = T
                    ) {
  cat('Saving ENM results \n')
  # Directories creation
  dir.create(path = paste0(directory, "/Rasters"))
  dir.create(path = paste0(directory, "/Tables"))
  
  # Raster saving
  cat('   rasters ...')
  setwd(paste0(directory, "/Rasters"))
  writeRaster(enm@proj[[1]], 'Probability', 'GTiff', overwrite = T)
  writeRaster(enm@proj[[2]], 'Niche', 'GTiff', overwrite = T)
  writeRaster(enm@incert, 'Uncertainity', 'GTiff', overwrite = T)
  cat('saved \n')
  
  # Tables saving
  cat('   tables ...')
  setwd(paste0(directory, "/Tables"))
  write.csv(enm@eval, 'ENMeval')
  write.csv(enm@algo.eval, 'AlgoEval')
  write.csv(enm@algo.corr, 'AlgoCorr')
  write.csv(enm@axes.contrib, 'VarImp')
  setwd(directory)
  cat('saved \n \n')
}

# 9 - functions forecasting #
ENM.modelling = function(Occurences, Env,
                         models = 'all',
                         # Data format
                         Spcol = 'SpeciesID', 
                         Xcol = 'Longitude',
                         Ycol = 'Latitude',
                         Norm = T,
                         # Modelling parameters
                         PA = 'Adapt', 
                         test = 'AIC', 
                         epsilon = 1e-08, 
                         maxit = 500, 
                         cv = 3, 
                         trees = 2500, 
                         thresh.shrink = 1e-03, 
                         final.leave = 1, 
                         split = 70,
                         # Projection parameters
                         AUC.thresh = 0.75, 
                         thresh = 1001,
                         # Saving
                         save = F,
                         directory = getwd()
                         ) {
  Occurences = TreatOcc(Occurences, Spcol, Xcol, Ycol)
  Env = TreatVar(Env, Norm = Norm)
  Occurences = Resample(Occurences, Env)
  Models = Modelling(Occurences, Env, models, PA, test, epsilon, maxit, cv, trees, thresh.shrink, final.leave, split)
  if (length(Models) > 0) {Models = Project(Models)}
  if (length(Models) > 0) {Models = AlgoENM(Models, AUC.thresh, thresh)}
  if (length(Models) > 0) {enm = ENMproj(Models, thresh)}
  if (save && exists('enm')) {save.enm(enm, directory)}
  if(exists('enm')) {
    return(enm)
  } else {
    cat('Ensemble modeling failed \n\n')
    return(NULL)
  }
}

# 10 - Species loop
sp.loop = function(Occurences, Env, 
                   models = c('all'),
                   save = T, 
                   name = format(Sys.time(), "%d.%m_%H:%M:%S"), 
                   log = T,
                   ...) {
  
  # Log file creation
  if (log) {sink(paste0('log_',format(Sys.time(), "%d.%m_%H:%M:%S")), append = T, type = c('output','message'), split = T)} 
  
  # Initializing log
  cat('################################################################################################################ \n')
  cat('New modeling (', as.character(models),') at',format(Sys.time(), "%H:%M:%S"),'\n')
  cat('################################################################################################################ \n\n')
  
  # Initianalizing saving
  directory = getwd()
  if (save) {dir.create(path = paste0(directory, "/",'Results_', name))}
  enm.list = list()
  
  # Launching modelling
  for (i in 1:length(levels(Occurences$SpeciesID))) {
    enm = ENM()
    Occurences_sp = subset(Occurences, Occurences$SpeciesID == levels(Occurences$SpeciesID)[i])
    if (save) {dir.create(path = paste0(directory, "/",'Results_', name, "/",levels(Occurences$SpeciesID)[i],'_Results'))}
    # Introduction de snowfall for a faster computing when all algorithms are ok
    enm = ENM.modelling(Occurences_sp, Env, models = models, ..., save = save, directory = paste0(directory, "/",'Results_', name, "/",levels(Occurences$SpeciesID)[i],'_Results'))
    enm.list[i] = enm
    
    # Saving ensemble model time in log
    cat('################################################################################################################ \n')
    cat(levels(Occurences$SpeciesID)[i],'done at',format(Sys.time(), "%H:%M:%S"),'\n')
    cat('################################################################################################################ \n\n')
    setwd(directory)
  }
  
  # Saving total time in log
  cat('################################################################################################################ \n')
  cat('Modeling done at',format(Sys.time(), "%H:%M:%S"),'\n')
  cat('################################################################################################################ \n\n')
  
  unlink('ID.*', recursive = T)
  if(log) {sink()}
  return(enm.list)
}

