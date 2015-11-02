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
#Occurences = subset(Occurences, Occurences$SpeciesID == c('6313','6314','6318'))

#### Main #### ----
sink('Log', split = T)
Sys.time()
# droplevels(Occurences[1:61,])
stack = Stack.Modelling('all', Occurences,  Env, rep = 1, tmp = F, save = F,
                        ensemble.metric = c('AUC','Kappa'), ensemble.thresh = c(0.7,0.2))
Sys.time()
sink()
plot(stack)

### Test all ### ----
GLM = Modelling('GLM',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
GAM = Modelling('ANN',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
enm1 = ensemble(GLM, GAM, uncertainity = F, name = '6318',
                ensemble.metric = c('AUC','Kappa','specificity'), ensemble.thresh = c(0.2,0,0.5), weight = F)
enm2 = Ensemble.Modelling(c('RF','ANN'), subset(Occurences, Occurences$SpeciesID == '6320'), Env, uncertainity = F)
stack1 = stacking(enm1, enm2)
plot(stack1)
stack2 = Stack.Modelling(c('GLM','GAM'), subset(Occurences, Occurences$SpeciesID == c('6318','6339')), Env, rep = 1, 
                         ensemble.metric = c('AUC','Kappa'), ensemble.thresh = c(0,0), weight = T, axes.metric = 'specificity')
plot(stack2)
