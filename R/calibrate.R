
# Calibration function based on Naimi's calibration function from the sdm-package (2016, version: 1.0-46)
# returns matrix with two columns, one for predicted probability & one for the proportion of presence per bin. The number of bins can be set through nbin.
# If plot==T, the result will be plotted as a calibration curve for quick visualization of the calibration success

calibrate <- function(x,nbin,plot){
  d <- data.frame(matrix(NA,nrow=nbin,ncol=2))
  colnames(d) <- c('predicted_probability','proportion_of_presences')
  p <-extract(x@projection,x@data[,c(1,2)])
  o <- x@data$Presence
  r <- max(p) - min(p)
  brk <- r / nbin
  for (i in 1:nbin){
    b <- min(p) + (i - 1) * brk 
    if(i < nbin) w <- which(p >= b & p < (b+brk))
    else w <- which(p >= b & p <= (b+brk))
    
    if(length(w) > 0) {
      d[i,2] <-  length(which(o[w] == 1)) / length(o[w])
      d[i,1] <- b+(brk/2)
    }
  }
  if(plot==T){
    plot(NA,NA, main=x@name,xlab="Predicted Probability of Occurrence [%]",ylab="Proportion of Presences [%]",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
    axis(2, at=seq(0,1,0.1), lab=paste0(seq(0,1,0.1) * 100), las=3)
    axis(1, at=seq(0,1,0.1), lab=paste0(seq(0,1,0.1) * 100))
    lines(d$predicted_probability,d$proportion_of_presences,col="coral")
    lines(seq(0,1,0.1),seq(0,1,0.1),lty=2)
  }
  return(d)
}

