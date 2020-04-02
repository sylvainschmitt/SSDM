# function SDMTools::accuracy() and related functions (auc, omission,
# sensitivity, specificity, prop.correct, Kappa)  removed from CRAn as not
# maintained

.accuracy <- function(obs,pred,threshold=0.5){
  #input checks
  if (length(obs)!=length(pred)) stop('this requires the same number of observed & predicted values')

  #deal with NAs
  if (length(which(is.na(c(obs,pred))))>0) {
    na = union(which(is.na(obs)),which(is.na(pred)))
    warning(length(na),' data points removed due to missing data')
    obs = obs[-na]; pred = pred[-na]
  }

  #define the n's and do checks
  n = length(obs); if (length(which(obs %in% c(0,1)))!=n) stop('observed values must be 0 or 1') #ensure observed are values 0 or 1

  # check / setup the threshold values
  if (length(threshold)==1 & threshold[1]<=1 & threshold[1]>=0) {
    thresholds = threshold
  } else if (length(threshold)==1 & threshold[1]>1) {
    thresholds = seq(0,1,length=threshold)
  } else if (length(threshold)>1 & max(threshold)<=1 & min(threshold)>=0) {
    thresholds = threshold
  } else { stop('inappropriate threshold values used as input. See help file.') }

  #cycle through each of the helpfiles
  out = data.frame(threshold=as.double(thresholds),AUC=NA,omission.rate=NA,sensitivity=NA,
                   specificity=NA,prop.correct=NA,Kappa=NA)
  for (ii in 1:length(thresholds)) {
    threshold = thresholds[ii]
    #convert preditions to binary based on threshold
    bin.pred = pred;
    if (threshold==0) {
      bin.pred[which(bin.pred>threshold)] = 1; bin.pred[which(bin.pred<=threshold)] = 0
    } else {
      bin.pred[which(bin.pred>=threshold)] = 1; bin.pred[which(bin.pred<threshold)] = 0
    }
    #create the confusion matrix ...just like using confusion.matrix command
    mat = table(bin.pred=factor(bin.pred,levels=c(0,1)),obs=factor(obs,levels=c(0,1)))
    attr(mat,'class') = 'confusion.matrix'
    #sppend data to output dataframe
    out$AUC[ii] = .auc(obs,bin.pred)
    out$omission.rate[ii] = .omission(mat)
    out$sensitivity[ii] = .sensitivity(mat)
    out$specificity[ii] = .specificity(mat)
    out$prop.correct[ii] = .prop.correct(mat)
    out$Kappa[ii] = .Kappa(mat)
  }
  #return the output data
  return(out)
}

.auc <- function(obs,pred){
  #input checks
  if (length(obs)!=length(pred)) stop('this requires the same number of observed & predicted values')


  #deal with NAs
  if (length(which(is.na(c(obs,pred))))>0) {
    na = union(which(is.na(obs)),which(is.na(pred)))
    warning(length(na),' data points removed due to missing data')
    obs = obs[-na]; pred = pred[-na]
  }

  #define the n's and do checks
  n = length(obs); if (length(which(obs %in% c(0,1)))!=n) stop('observed values must be 0 or 1') #ensure observed are values 0 or 1
  n1 = as.double(length(which(obs==1))); n0 = as.double(length(which(obs==0)))
  if (n1==0 || n1==n) return( NaN ) #if all observed 1's or 0's return NaN

  ###calc AUC
  pred0 = pred[which(obs==0)]
  pred1 = pred[which(obs==1)]
  ranks = rank(pred,ties.method='average')#define ranks
  ranks0 = ranks[which(obs==0)]
  ranks1 = ranks[which(obs==1)]
  U = n0*n1 + (n0*(n0+1))/2 - sum(ranks0) #calc U stat
  AUC = U/(n0*n1) #estimate AUC
  if (AUC<.5) AUC = 1-AUC

  #return the auc value
  return(AUC)
}

.omission <- function(mat){
  #input checks
  if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
  #return the value
  return(mat[1,2]/sum(mat[,2]))
}

.sensitivity = function(mat) {
  #input checks
  if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
  #return the value
  return(mat[2,2]/sum(mat[,2]))
}

.specificity = function(mat) {
  #input checks
  if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
  #return the value
  return(mat[1,1]/sum(mat[,1]))
}

.prop.correct = function(mat) {
  #input checks
  if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
  #return the value
  return(sum(diag(mat))/sum(mat))
}

.Kappa <- function(mat){
  #input checks
  if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
  #calculate Kappa
  n<-sum(mat)
  colsums<-as.double(apply(mat,2,sum))
  rowsums<-as.double(apply(mat,1,sum))
  t1 = sum(diag(mat))/n; 	t2 = sum(rowsums*colsums)/(n^2)
  #return the value
  return((t1-t2)/(1-t2))
}
