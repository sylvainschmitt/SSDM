#### Loading #### ----
setwd("/home/amap/Documents/R/Sans Biomod 2/")
source('Data.prep.R')
source('Classes.R')
source('Modelling.R')
Env = load.var('./data', format = c('tif'), factors = c('substrat', 'secteurs'))
Env = stack(aggregate(Env, 10, fun = max))
Env[[7]] = as.factor(Env[[7]])
Env[[8]] = as.factor(Env[[8]])
Occurences = load.occ('data', Env, Spcol = 'SpeciesID', dec = '.')
#Occurences = subset(Occurences, Occurences$SpeciesID == c('6339','6313','6318'))

#### Main #### ----
sink('Log', split = T)
Sys.time()
stack = Stack.Modelling('all', Occurences, Env, rep = 10, tmp = T, save = T)
Sys.time()
sink()
plot(stack)

### Test all ### ----
GLM = Modelling('GLM',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
GAM = Modelling('GAM',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
enm1 = ensemble(GLM, GAM, uncertainity = F, AUCthresh = 0)
enm2 = Ensemble.Modelling(c('RF','ANN'), subset(Occurences, Occurences$SpeciesID == '6320'), Env, uncertainity = F, AUCthresh = 0)
stack1 = stacking(enm1, enm2)
stack2 = Stack.Modelling(c('SVM','ANN'), subset(Occurences, Occurences$SpeciesID == c('6318','6319')), Env, rep = 3)
plot(stack2)
