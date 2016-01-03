.checkargs = function(Xcol = 'Longitude',
                      Ycol = 'Latitude',
                      Pcol = NULL,
                      Spcol = 'SpeciesID',
                      name = NULL,
                      Spname = NULL,
                      save = F,
                      path = getwd(),
                      PA = NULL,
                      rep = 1,
                      cv = 'holdout',
                      cv.param = c(0.7,2),
                      thresh = 1001,
                      metric = 'SES',
                      axes.metric = 'Pearson',
                      select = F,
                      select.metric = c('AUC'),
                      select.thresh = c(0.75),
                      verbose = T,
                      GUI = F,
                      uncertainty = T,
                      tmp = F,
                      ensemble.metric = c('AUC'),
                      ensemble.thresh = c(0.75),
                      weight = T,
                      method = 'P',
                      rep.B = 1000,
                      GeoRes = T,
                      reso = 1,
                      file = NULL,
                      files = NULL,
                      format = c('.grd','.tif','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img'),
                      categorical = NULL,
                      Norm = T,
                      enm = Ensemble.SDM(),
                      stack = Stacked.SDM(),
                      range = NULL,
                      endemism = 'WEI',
                      cores = 1){
  ## Argument checking function ##
  # Occurrences, Environment, and part of X, Y, Pcol are directly defined in functions

  # Col
  if(!inherits(Xcol, 'character')){stop('Xcol parameter should be a character (column name).')}
  if(!inherits(Ycol, 'character')){stop('Ycol parameter should be a character (column name).')}
  if(!inherits(Pcol, 'character') && !is.null(Pcol)){stop('Pcol parameter should be a character (column name), or NULL if no presence column.')}
  if(!inherits(Spcol, 'character') && !is.null(Pcol)){stop('Spcol parameter should be a character (column name), or NULL if no species column..')}

  # Name and Spname
  if(!inherits(name, 'character') && !is.null(name)){stop('name parameter should be a character, or NULL to be automatically defined.')}
  if(!inherits(Spname, 'character') && !is.null(Spname)){stop('Spname parameter should be a character, or NULL to be automatically defined.')}

  # PA
  if(!inherits(PA,'list') && !is.null(PA)){stop('PA parameter should be a list, or NULL to be automatically defined.')}
  if(!is.null(PA) && (is.null(PA$nb) || is.null(PA$strat))) {stop('PA should be a list containing as PA = list(\'nb\'=..., \'strat\'=...).')}
  if(!is.null(PA) && !is.null(PA$nb) && (abs(PA$nb-round(PA$nb))!=0 || PA$nb < 1)){stop('PA$nb should be an integer > 0.')}
  if(!is.null(PA) && !is.null(PA$strat) && !(PA$strat %in% c('random','disk'))){stop('PA$strat should be random or disk (see help).')}

  # rep
  if(!inherits(rep,'numeric') || abs(thresh -round(thresh))!=0 || thresh < 1){stop('rep parameter should be an integer > 1.')}

  # CV
  if(!inherits(cv, 'character')){stop('cv parameter should be a character.')}
  if(!(cv %in% c('k-fold','holdout','LOO'))){stop('cv parameter should be k-fold, holdout or LOO (see help).')}
  if(!inherits(cv.param, 'numeric') || is.null(cv.param)){stop('cv.param parameter should be a numeric')}
  if(cv != 'LOO' && length(cv.param) < 2){stop('cv.param parameter should be of length 2')}
  if(cv == 'k-fold' && ((abs(cv.param[1]-round(cv.param[1]))!=0 || abs(cv.param[2]-round(cv.param[2]))!=0) || (cv.param[1]==0 || cv.param[2]==0))){stop('cv.param parameters (k and repetitions) should be both integers > 0 (see help).')}
  if(cv == 'holdout' && (abs(cv.param[2]-round(cv.param[2]))!=0 || cv.param[2]==0)){stop('cv.param[2] (repetitions) parameters should be an integer > 0 (see help).')}
  if(cv == 'holdout' && (cv.param[1]<0 || cv.param[1]>1)){stop('cv.param[1] (train fraction) parameters should be a float between 0 and 1 (see help).')}

  # thresh
  if(!inherits(thresh,'numeric') || abs(thresh -round(thresh))!=0){stop('thresh parameter should be an integer.')}

  # metric
  if(!inherits(metric,'character') || !(metric %in% c('Kappa','CCR','TSS','SES','LW','ROC'))){stop('metric parameter should be Kappa, CCR, TSS, SES, LW, or ROC (see help).')}

  # axes.metric
  if(!inherits(axes.metric,'character') || !(axes.metric %in% c('Pearson','AUC','Kappa','omission.rate','sensitivity','specificity','prop.correct'))){stop('axes.metric parameter should be Pearson, AUC, Kappa, omission.rate, sensitivity, specificity or prop.correct (see help).')}

  # select
  if(!inherits(select,'logical')){stop('select parameter should be a logical (True or False).')}
  if(!inherits(select.metric,'character')){stop('select.metric parameter should be a character.')} else {
    for (i in 1:length(select.metric)){
      if(!(select.metric[i] %in% c('AUC','Kappa','sensitivity','specificity','prop.correct'))){stop(paste('select.metric',i,'parameter should be AUC, Kappa, sensitivity, specificity, or prop.correct (see help).'))}
    }
  }
  if(!inherits(select.thresh,'numeric')){stop('select.thresh parameter should be a numeric.')} else {
    for (i in 1:length(select.thresh)){
      if(select.thresh[i]<0 || select.thresh[i]>1){stop(paste('select.thresh',i,'parameter should be a float between 0 and 1.'))}
    }
  }
  if(length(select.thresh) != length(select.metric)){stop('select.thresh and select.metric parameters should have the same length to correspond (see help).')}

  # Verbose, GUI, uncertainty, tmp and save
  if(!inherits(verbose,'logical')){stop('verbose parameter should be a logical (True or False).')}
  if(!inherits(GUI,'logical')){stop('GUI parameter should be a logical (True or False).')}
  if(!inherits(uncertainty,'logical')){stop('uncertainty parameter should be a logical (True or False).')}
  if(!inherits(tmp,'logical')){stop('tmp parameter should be a logical (True or False).')}
  if(!inherits(save,'logical')){stop('save parameter should be a logical (True or False).')}
  if(!inherits(path, 'character')  && !is.null(path)){stop('path parameter should be a character.')}

  # Ensemble
  if(!inherits(weight,'logical')){stop('weight parameter should be a logical (True or False).')}
  if(!inherits(ensemble.metric,'character')){stop('ensemble.metric parameter should be a character.')} else {
    for (i in 1:length(ensemble.metric)){
      if(!(ensemble.metric[i] %in% c('AUC','Kappa','sensitivity','specificity','prop.correct'))){stop(paste('ensemble.metric',i,'parameter should be AUC, Kappa, sensitivity, specificity, or prop.correct (see help).'))}
    }
  }
  if(!inherits(ensemble.thresh,'numeric')){stop('ensemble.thresh parameter should be a numeric.')} else {
    for (i in 1:length(ensemble.thresh)){
      if(ensemble.thresh[i]<0 || ensemble.thresh[i]>1){stop(paste('ensemble.thresh',i,'parameter should be a float between 0 and 1.'))}
    }
  }
  if(length(ensemble.thresh) != length(ensemble.metric)){stop('select.thresh and select.metric parameters should have the same length to correspond (see help).')}

  # Stacking
  if(!inherits(method,'character') || !(method %in% c('P','B','T'))){stop('method parameter should be P, B, or T (see help).')}
  if(method == 'B' && !inherits(rep.B,'numeric') && abs(rep.B-round(rep.B)) != 0 && rep.B < 1){stop('rep.B parameter should be an integer > 1 (see help).')}

  # load.occ and load.var
  if(!inherits(GeoRes,'logical')){stop('GeoRes parameter should be a logical (True or False).')}
  if(!inherits(reso,'numeric')){stop('reso parameter should be numeric.')}
  if(!inherits(file,'character') && !is.null(file)){stop('file parameter should be characters or NULL')}
  if(!inherits(files,'character') && !is.null(files)){stop('files parameter should be characters or NULL')}
  if(!is.null(format) && !inherits(format,'character') || (inherits(format,'character') && !(format %in% c('.grd','.tif','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img')))){stop('format parameter should be .grd, .tif, .asc, .sdat, .rst, .nc, .tif, .envi, .bil or .img')}
  if(!inherits(categorical,'character') && !is.null(categorical)){stop('categorical parameter should be characters or NULL')}
  if(!inherits(Norm,'logical')){stop('Norm parameter should be a logical (True or False).')}

  # save
  if(!inherits(enm,'Ensemble.SDM')){stop('enm parameter should be an Ensemble.SDM.')}
  if(!inherits(stack,'Stacked.SDM')){stop('stack parameter should be a Stacked.SDM.')}

  # Range and Endemism parameters
  if(!is.null(endemism)){if(!inherits(endemism,'character') || !(endemism %in% c('Any','WEI','CWEI'))){stop('endemism parameter should be Any, WEI, or CWEI (see help).')}}
  if(!is.null(range)){
    if(!inherits(range,'numeric') && range < 0){
      stop('range parameter should be numeric and > 0 (see help).')
    }
  }

  # Cores
  if(!inherits(cores,'numeric') || abs(cores-round(cores)) != 0 || cores < 1){stop('cores parameter should be an integer > 0 (see help).')}
}
