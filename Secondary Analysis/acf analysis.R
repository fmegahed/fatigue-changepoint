rm(list = ls()) # clear global environment
graphics.off() # close all graphics
library(pacman) # needs to be installed first
p_load(R.matlab, plotly, extrafont, grDevices, gridExtra,
       dplyr, stringr, tidyverse, utils, reshape2,
       anomalize, MVN, forecast, fractal,
       ecp, dfphase1, 
       MALDIquant, TSclust,
       knitr, kableExtra, dplyr)

load(file = "./Data/RGenerated/ScaledData.RData")

med.fil<-function(x,d,w){
  ms = {}
  i = 1
  while( w+d*(i-1) < length(x))
  {
    ms[i] = median(x[(1+d*(i-1)):(w+d*(i-1))])
    i = i + 1
  }
  return(ms)
}

percent.time.from.start <- med.fil(subject10_features_norm[,1], d=21, w=21)
scaled.stride.len <- med.fil(subject10_features_norm[,2], d=21, w=21)
# scaled.stride.height <- med.fil(subject10_features_norm[,3], d=21, w=21)
# stride.duration <- med.fil(subject10_features_norm[,4], d=21, w=21)
sub10.nooverlap <- data.frame(cbind(percent.time.from.start, scaled.stride.len))

# sub10.nooverlap <- data.frame(cbind(percent.time.from.start, scaled.stride.len, 
#                                     scaled.stride.height, stride.duration))

percent.time.from.start <- med.fil(subject10_features_norm[,1], d=11, w=21)
scaled.stride.len <- med.fil(subject10_features_norm[,2], d=11, w=21)
sub10.Halfoverlap <- data.frame(cbind(percent.time.from.start, scaled.stride.len))

percent.time.from.start <- med.fil(subject10_features_norm[,1], d=1, w=21)
scaled.stride.len <- med.fil(subject10_features_norm[,2], d=1, w=21)
sub10.wMinus1overlap <- data.frame(cbind(percent.time.from.start, scaled.stride.len))


z1= acf(sub10.nooverlap$scaled.stride.len, xlab ="Lag", ylab="ACF", main="ACF P10's Stride Length (no overlap)")

z2= acf(sub10.Halfoverlap$scaled.stride.len, xlab ="Lag", ylab="ACF", main="ACF Plot for P10's Stride Length (overlap = w/2)")

z3= acf(sub10.wMinus1overlap$scaled.stride.len, xlab ="Lag", ylab="ACF", main="ACF Plot for P10's Stride Length (overlap = w-1)")

p1 <- plot(sub10.nooverlap$percent.time.from.start, sub10.nooverlap$scaled.stride.len,
           main="No overlap Median Filter for P10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l") 

p2 <- plot(sub10.Halfoverlap$percent.time.from.start, sub10.Halfoverlap$scaled.stride.len,
           main="W/2 overlap Median Filter for P10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l")

p3 <- plot(sub10.wMinus1overlap$percent.time.from.start, sub10.wMinus1overlap$scaled.stride.len,
           main="W/2 overlap Median Filter for P10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l")

windows()
par(mfrow= c(2, 3))

plot(sub10.nooverlap$percent.time.from.start, sub10.nooverlap$scaled.stride.len,
           main="No overlap Median Filter for Participant 10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l", cex.lab=1.5, cex=1.5) 

plot(sub10.Halfoverlap$percent.time.from.start, sub10.Halfoverlap$scaled.stride.len,
           main="W/2 overlap Median Filter for Participant 10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l", cex.lab=1.5, cex=1.5)

plot(sub10.wMinus1overlap$percent.time.from.start, sub10.wMinus1overlap$scaled.stride.len,
           main="W/2 overlap Median Filter for Participant 10", xlab="Percent time from start",
           ylab="Scaled Stride Length", type="l", cex.lab=1.5, cex=1.5)

acf(sub10.nooverlap$scaled.stride.len, xlab ="Lag", ylab="ACF", main="ACF P10's Stride Length (no overlap)", 
    cex.lab=1.5, cex=1.5, lag.max = 30)

acf(sub10.Halfoverlap$scaled.stride.len, xlab ="Lag", ylab="ACF", 
    main="ACF Plot for P10's Stride Length (overlap = w/2)", cex.lab=1.5, 
    cex=1.5, lag.max = 30)

acf(sub10.wMinus1overlap$scaled.stride.len, xlab ="Lag", ylab="ACF", 
    main="ACF Plot for P10's Stride Length (overlap = w-1)", cex.lab=1.5, 
    cex=1.5, lag.max = 30)

