scratch <- read.csv("scratch.txt", head=FALSE)
n <- nrow(scratch)
scratch2 <- data.frame(t(rep("", 10)), stringsAsFactors=FALSE)

for(i in 1:n){
  splits <- strsplit(as.character(scratch[i,]), " ")[[1]]
  k <- length(splits)
  if(k == 9)
      splits <- c(splits, "")
  scratch2[i,] <- splits
}

dat <- scratch2[,c(3,5,7,9,10)]
colnames(dat) <- c("sampler", "T", "V", "W", "Finish")
datstart <- dat[dat$Finish=="",]
datfin <- dat[dat$Finish=="FINISHED",]

datstart$found <- 0
datfin$found <- 0

ks <- nrow(datstart)
kf <- nrow(datfin)

for(i in 1:ks){
  for(j in 1:kf){
    if(datfin$sampler[j]==datstart$sampler[i]){
      if(datfin$T[j]==datstart$T[i]){
        if(datfin$V[j]==datstart$V[i]){
          if(datfin$W[j]==datstart$W[i]){
            datfin$found[j] <- 1
            datstart$found[i] <- 1
          }
        }
      }
    }
  }
}

datstart[datstart$found==0,]
datfin[datfin$found==0,]

sampler    T      V       W Finish found
596   sdalt T=10 QV=100 QW=1000            0

   sampler      T       V       W   Finish found
7    sdalt  T=100 QV=0.01 QW=0.01 FINISHED     0
19   sdalt T=1000 QV=0.01 QW=0.01 FINISHED     0

for(i in 1:ks){
  datstart$T[i] <- strsplit(datstart$T[i], "T=")[[1]][2]
  datstart$V[i] <- strsplit(datstart$V[i], "QV=")[[1]][2]
  datstart$W[i] <- strsplit(datstart$W[i], "QW=")[[1]][2]
}
for(i in 1:kf){
  datfin$T[i] <- strsplit(datfin$T[i], "T=")[[1]][2]
  datfin$V[i] <- strsplit(datfin$V[i], "QV=")[[1]][2]
  datfin$W[i] <- strsplit(datfin$W[i], "QW=")[[1]][2]
}

datfin$T <- as.numeric(datfin$T)
datfin$V <- round(as.numeric(datfin$V),4)
datfin$W <- round(as.numeric(datfin$W),4)
datstart$T <- as.numeric(datstart$T)
datstart$V <- round(as.numeric(datstart$V),4)
datstart$W <- round(as.numeric(datstart$W),4)

datfin10 <- table(datfin[datfin$T==10,c(3,4)])
datfin100 <- table(datfin[datfin$T==100,c(3,4)])
datfin1000 <- table(datfin[datfin$T==1000,c(3,4)])
datstart10 <- table(datstart[datstart$T==10,c(3,4)])
datstart100 <- table(datstart[datstart$T==100,c(3,4)])
datstart1000 <- table(datstart[datstart$T==1000,c(3,4)])


idf10 <- which(datfin10==0)
idf100 <- which(datfin100==0)
idf1000 <- which(datfin1000==0)
ids10 <- which(datstart10==0)
ids100 <- which(datstart100==0)
ids1000 <- which(datstart1000==0)

rf10 <- idf10%%11
rf10[rf10==0] <- 11
cf10 <- (idf10+10)%/%11
rf100 <- idf100%%11
rf100[rf100==0] <- 11
cf100 <- (idf100+10)%/%11
rf1000 <- idf1000%%11
rf1000[rf1000==0] <- 11
cf1000 <- (idf1000+10)%/%11

rs10 <- ids10%%11
rs10[rs10==0] <- 11
cs10 <- (ids10+10)%/%11
rs100 <- ids100%%11
rs100[rs100==0] <- 11
cs100 <- (ids100+10)%/%11
rs1000 <- ids1000%%11
rs1000[rs1000==0] <- 11
cs1000 <- (ids1000+10)%/%11

V <- 10^(c(0:10)/2-2)
W <- V

s10m <- data.frame(V=V[rs10], W=W[cs10], SF="S", T=10)
s100m <- data.frame(V=V[rs100], W=W[cs100], SF="S", T=100)
s1000m <- data.frame(V=V[rs1000], W=W[cs1000], SF="S", T=1000)
f10m <- data.frame(V=V[rf10], W=W[cf10], SF="F", T=10)
f100m <- data.frame(V=V[rf100], W=W[cf100], SF="F", T=100)
f1000m <- data.frame(V=V[rf1000], W=W[cf1000], SF="F", T=1000)

rbind(s10m, s100m, s1000m, f10m, f100m, f1000m)

V[unique(c(rs10, rs100, rs1000, rf10, rf100, rf1000))]
W[unique(c(cs10, cs100, cs1000, cf10, cf100, cf1000))]

