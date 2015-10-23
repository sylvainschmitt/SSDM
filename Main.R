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

stack = Stack.Modelling(c('GLM', 'RF'), Occurences, Env, tmp = T, AUCthresh = 0, save = T)
plot(stack)
