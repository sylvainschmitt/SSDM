#### Loading #### ----
setwd("/home/amap/Documents/R/Sans Biomod 2/")
source('Data.prep.R')
source('Classes.R')
source('Modelling.R')
Env = load.var('./data', format = 'tif', factors = c('substrat', 'secteurs'))
Env = stack(aggregate(Env, 10, fun = max))
Env[[7]] = as.factor(Env[[7]])
Env[[8]] = as.factor(Env[[8]])
Occurences = load.occ('data', Env, dec = '.')
Occurences = subset(Occurences, Occurences$SpeciesID == '6339')

#### Main #### ----
# enm = Ensemble.Modelling(c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'),
#                          Occurences, Env, name = '6339', rep = 10, save = T, tmp = T)

enm1 = Ensemble.Modelling(c('GLM','GAM','MARS'), Occurences, Env, name = '6339', AUCthresh = 0)
enm2 = Ensemble.Modelling(c('GLM','MARS','GAM'), Occurences, Env, name = '6339', AUCthresh = 0)
stack = stacking(enm1, enm2)
plot(stack)
