source("varfun.R")

dat2010 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2010.csv")
dat2011 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2011.csv")

## finds the 3 heated 25cm temperature 2011 time series and averages them
idx <- grep("H.25cm.T", colnames(dat2011))
datheat2011 <- dat2011[,idx]
datMH11 <- apply(datheat2011,1,mean)

## Fit model 3b using MLE
## note: parameter estimates are on log scale
MLEest11 <- dlmMLE(datMH11, rep(0,4), mymodMLE, other=c(2,24,3))
MLEstart <- exp(MLEest11$par) # my estimates don't replicate


