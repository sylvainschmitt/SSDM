### Loading #### ----
directory = '/home/amap/Documents/data'
Env = load.var(directory, format = c('tif'), factors = c('SUB', 'PHY'))
library(raster)
Env = stack(aggregate(Env, 10, fun = max))
Env[[3]] = as.factor(Env[[3]])
Env[[6]] = as.factor(Env[[6]])
Occurences = load.occ(directory, Env, Spcol = 'SpeciesID', dec = '.')
#Occurences = subset(Occurences, Occurences$SpeciesID == c('6313','6314','6318'))

#### Main #### ----
sink('Log', split = T)
Sys.time()
# droplevels(Occurences[1:61,])
stack = Stack.Modelling('all', Occurences,  Env, rep = 10, tmp = F, save = T,
                        ensemble.metric = c('AUC','Kappa'), ensemble.thresh = c(0.7,0.2))
Sys.time()
sink()
plot(stack)

### Test all ### ----
GLM = Modelling('GLM',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
GAM = Modelling('RF',subset(Occurences, Occurences$SpeciesID == '6318'), Env)
enm1 = ensemble(GLM, GAM, uncertainity = F, name = '6318', ensemble.metric = c('AUC','Kappa'), ensemble.thresh = c(0,0))
enm2 = Ensemble.Modelling(c('RF','ANN'), subset(Occurences, Occurences$SpeciesID == '6319'), Env, uncertainity = F)
stack1 = stacking(enm1, enm2)
plot(stack1)
stack2 = Stack.Modelling(c('ANN','RF'), droplevels(Occurences[1:61,]), Env, rep = 1, ensemble.metric = c('AUC','Kappa'),
                         ensemble.thresh = c(0.6,0.1), weight = T, axes.metric = 'specificity')
plot(stack2)
