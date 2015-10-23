#### Loading #### ----
setwd("/home/amap/Documents/R/Sans Biomod 2/")
source('Data.prep.R')
source('Classes.R')
source('Modelling.R')
Env = load.var('./data', format = c('tif'), factors = c('substrat', 'secteurs'))
Env = stack(aggregate(Env, 10, fun = max))
Env[[7]] = as.factor(Env[[7]])
Env[[8]] = as.factor(Env[[8]])
Occurences = load.occ('data', Env, dec = '.')
Occurences = subset(Occurences, Occurences$SpeciesID == c('6339','6313','6318'))

#### Main #### ----
# stack = Stack.Modelling(c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'), Occurences, Env, rep = 10, tmp = T)

enm1 = Ensemble.Modelling(c('GLM', 'RF'), subset(Occurences, Occurences$SpeciesID == '6318'), 
                          name = '6318', Env, tmp = T, AUCthresh = 0)
enm2 = Ensemble.Modelling(c('GLM', 'RF'), subset(Occurences, Occurences$SpeciesID == '6313'), 
                          name = '6313', Env, tmp = T, AUCthresh = 0)

stack = stacking(enm1, enm2)
plot(stack)
