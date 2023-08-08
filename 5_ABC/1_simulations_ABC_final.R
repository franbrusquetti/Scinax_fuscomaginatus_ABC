## in order to use this code you need to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

## set working directory
setwd(“home/Documents/Simulation_Data”)

### variable declarations

## number of simulations
numsim <- 500000
# mitochondrial
Nsammt = 60
## sample size of Pop1
PopA = 35
## sample size of Pop2
PopB = 20
## sample size of Pop3
PopC = 3
## sample size of Pop4
PopD = 2
## number of base pairs
Lmt <- 615
## number of years per generation
genlen <- 1
## mean mutation rate for COI-amphibians
mutratemt <- 4.305*10^(-6)
# bfib
Nsam1 = 62*2
PopA1 = 36*2
PopB1 = 21*2
PopC1 = 3*2
PopD1 = 2*2
# 27_28
Nsam2 = 57*2
PopA2 = 36*2
PopB2 = 18*2
PopC2 = 1*2
PopD2 = 2*2
# rpl3
Nsam3 = 40*2
PopA3 = 20*2
PopB3 = 15*2
PopC3 = 4*2
PopD3 = 1*2
## number of base pairs
Lnu <- 374 #average of nuclear genes
## number of loci (nuclear genes)
numloc <- 3

### define default values for all parameters that are specific to a subset of the models
# Set founded population size ratio to 1, equal to the ancestral population.
AncRatioPopA <- 1
AncRatioPopB <- 1
AncRatioPopC <- 1
AncRatioPopD <- 1
# Set current population size ratio to 1, equal to the ancestral population.
SizeRatioA <- 1
SizeRatioB <- 1
SizeRatioC <- 1
SizeRatioD <- 1
# Set growth rates to 0.
GrowthA <- 0
GrowthB <- 0
GrowthC <- 0
GrowthD <- 0
MigrationAB <- 0
MigrationAC <- 0
MigrationAD <- 0
MigrationBA <- 0
MigrationBC <- 0
MigrationBD <- 0
MigrationCA <- 0
MigrationCB <- 0
MigrationCD <- 0
MigrationDA <- 0
MigrationDB <- 0
MigrationDC <- 0

## create a variable to store all the parameters
parameters <- data.frame()

### vicariant model (model 1): constant population sizes after splits, with first split corresponding to Central Brazilian Plateau (CBP) uplift

for (i in 1:numsim) {
  
  ### Define parameters
  # based on parameters estimated under multinomial logistic regression (MNLOG) with abc function (tolerance 0.05), minimum and maximum values correspond to 95% Quantile interval
  Ne <- runif(1, 1313770.5703, 5976174.0839 )
  # based on parameters estimated under multinomial logistic regression (MNLOG) with abc function (tolerance 0.05), minimum and maximum values correspond to 95% Quantile interval
  thetamt <- runif(1, 1.9229, 121.1588)
  # based on parameters estimated under multinomial logistic regression (MNLOG) with abc function (tolerance 0.05), minimum and maximum values correspond to 95% Quantile interval
  thetanuc <- runif(1, 3.8535, 64.9567)
  # Prior for time of the most ancient event, follows a lognormal distribution with mean 14.9 and standard deviation 0.6,  which corresponds to mid-tertiary to late-tertiary/early-quaternary, 2 to 7 Mya, final plateau uplift (in years: 95% HPD= 1100000-7940000, median= 2960000)
  T3 <- rlnorm(1, 14.9, 0.6)
  # time in coalescent units
  coalT3 <- T3/Ne
  # Prior for time of an ancient event that is more recent than T3, follows a uniform distribution bounded 0.01 and T3.
  coalT2 <- runif(1, 0.01, coalT3)
  # Prior for time of an ancient event that is more recent than T2, follows a uniform distribution bounded 0.01 and T2.
  coalT1 <- runif(1, 0.01, coalT2)
  
  ## ms's command: 4 populations with stable demography, no migration
  # locus 1
  system(sprintf("./ms 124 1 -t %f -I 4 72 42 6 4 -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> vic.txt", thetanuc, coalT1, coalT2, coalT3))
  # locus 2
  system(sprintf("./ms 114 1 -t %f -I 4 72 36 2 4 -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> vic.txt", thetanuc, coalT1, coalT2, coalT3))
  # locus 3
  system(sprintf("./ms 80 1 -t %f -I 4 40 30 8 2 -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> vic.txt",  thetanuc, coalT1, coalT2, coalT3))
  # mitochondrial
  system(sprintf("./ms 60 1 -t %f -I 4 35 20 3 2 -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> vic_mit.txt",  thetamt, coalT1, coalT2, coalT3))
  
  
  ## save parameter values
  parameters <- rbind(parameters, data.frame(Ne, thetanuc, thetamt, coalT1, coalT2, coalT3))
}

## plateau-depression model (model 2), colonization scenario, with stable populations in plateaus and small ancient populations in depressions (founder effect), with first split corresponding to CBP uplift
for (i in 1:numsim) {
  
  ### Define parameters
  Ne <- runif(1, 1313770.5703, 5976174.0839 )
  thetamt <- runif(1, 1.9229, 121.1588)
  thetanuc <- runif(1, 3.8535, 64.9567)
  T3 <- rlnorm(1, 14.9, 0.6)
  coalT3 <- T3/Ne
  coalT2 <- runif(1, 0.01, coalT3)
  coalT1 <- runif(1, 0.01, coalT2)
  # Sample the ratio between the founded and the ancestral population sizes
  FoundedSizeRatioA <- runif(1, 0, 0.1)
  FoundedSizeRatioC <- runif(1, 0, 0.1)
  # Uses a beta distribution to sample the ratio between the current and the founded population sizes. The shape of the distribution results in a median value of 0.0452 and a 95% Quantile interval from 0.00169 to 0.218
  AncRatioPopA <- rbeta(1, 1.0, 15.0)
  AncRatioPopC <- rbeta(1, 1.0, 15.0)
  GrowthA = -(1/coalT2)*log(AncRatioPopA)
  GrowthC = -(1/coalT1)*log(AncRatioPopC)
  # Find the ratio between the current and the ancestral population sizes.
  SizeRatioA = FoundedSizeRatioA/AncRatioPopA
  SizeRatioC = FoundedSizeRatioC/AncRatioPopC
  
  ## ms's command: 4 populations with stable demography in populations B and D and founder effect in populations A and C, no migration
  # locus 1
  system(sprintf("./ms 124 1 -t %f -I 4 72 42 6 4 -g 3 %f -g 1 %f -en %f 3 %f -ej %f 3 1 -en %f 1 %f -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> p_d.txt", thetanuc, GrowthC, GrowthA, coalT1, FoundedSizeRatioC, coalT1, coalT2, FoundedSizeRatioA, coalT2, coalT3))
  # locus 2
  system(sprintf("./ms 114 1 -t %f -I 4 72 36 2 4 -g 3 %f -g 1 %f -en %f 3 %f -ej %f 3 1 -en %f 1 %f -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> p_d.txt", thetanuc, GrowthC, GrowthA, coalT1, FoundedSizeRatioC, coalT1, coalT2, FoundedSizeRatioA, coalT2, coalT3))
  # locus 3
  system(sprintf("./ms 80 1 -t %f -I 4 40 30 8 2 -g 3 %f -g 1 %f -en %f 3 %f -ej %f 3 1 -en %f 1 %f -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> p_d.txt", thetanuc, GrowthC, GrowthA, coalT1, FoundedSizeRatioC, coalT1, coalT2, FoundedSizeRatioA, coalT2, coalT3))
  # mitochondrial
  system(sprintf("./ms 60 1 -t %f -I 4 35 20 3 2 -g 3 %f -g 1 %f -en %f 3 %f -ej %f 3 1 -en %f 1 %f -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> p_d_mit.txt", thetamt, GrowthC, GrowthA, coalT1, FoundedSizeRatioC, coalT1, coalT2, FoundedSizeRatioA, coalT2, coalT3))
  
  ## save parameter values
  parameters <- rbind(parameters, data.frame(Ne, thetanuc, thetamt, coalT1, coalT2, coalT3))
}

## isolation by instability model (model 3), represents a refugia scenario with all populations experiencing size reductions, with first split corresponding to Last Interglacial (LIG)
for (i in 1:numsim) {
  
  ### Define parameters   
  Ne <- runif(1, 1313770.5703, 5976174.0839 )
  thetamt <- runif(1, 1.9229, 121.1588)
  thetanuc <- runif(1, 3.8535, 64.9567)
  # Prior for time of the most ancient event, follows a lognormal distribution with mean 12 and standard deviation 0.6,  which corresponds to (150000 years - LIG)(in years: 95% HPD= 60700-437000, median= 163000)
  T3 <-rlnorm(1, 12, 0.6)
  coalT3 <- T3/Ne
  coalT2 <- runif(1, 0.00, coalT3)
  coalT1 <- runif(1, 0.00, coalT2)
  ## bottleneck in all populations. Prior sampled from a beta distribution resulted in median= 0.264, 95% HPD from 0.0628 to 0.582
  SizeRatioA <- rbeta(1, 2.0, 5.0)
  SizeRatioB <- rbeta(1, 2.0, 5.0)
  SizeRatioC <- rbeta(1, 2.0, 5.0)
  SizeRatioD <- rbeta(1, 2.0, 5.0)
  
  ## ms's command: 4 populations with bottleneck in all populations, no migration
  # locus 1
  system(sprintf("./ms 124 1 -t %f -I 4 72 42 6 4 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ref1.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # locus 2
  system(sprintf("./ms 114 1 -t %f -I 4 72 36 2 4 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ref1.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # locus 3
  system(sprintf("./ms 80 1 -t %f -I 4 40 30 8 2 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ref1.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # mitochondrial
  system(sprintf("./ms 60 1 -t %f -I 4 35 20 3 2 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ref1_mit.txt", thetamt, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  
  ## save parameter values
  parameters <- rbind(parameters, data.frame(Ne, thetanuc, thetamt, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
}

## isolation by instability model (model 4), represents a refugia scenario but with stable populations in the CBP, with first split corresponding to LIG
for (i in 1:numsim) {
  
  #### Define parameters
  Ne <- runif(1, 1313770.5703, 5976174.0839 )
  thetamt <- runif(1, 1.9229, 121.1588)
  thetanuc <- runif(1, 3.8535, 64.9567)
  T3 <-rlnorm(1, 12, 0.6)
  coalT3 <- T3/Ne
  coalT2 <- runif(1, 0.00, coalT3)
  coalT1 <- runif(1, 0.00, coalT2)
  ## bottleneck in populations A and C. Prior sampled from a beta distribution resulted in median= 0.264, 95% HPD from 0.0628 to 0.582
  SizeRatioA <- rbeta(1, 2.0, 5.0)
  SizeRatioC <- rbeta(1, 2.0, 5.0)
  
  ## ms's command: 4 populations with bottleneck only in populations A and C, no migration
  # locus1
  system(sprintf("./ms 124 1 -t %f -I 4 72 42 6 4 -n 3 %f -n 1 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> ref2.txt", thetanuc, SizeRatioC, SizeRatioA, coalT1, coalT2, coalT3))
  # locus 2
  system(sprintf("./ms 114 1 -t %f -I 4 72 36 2 4 -n 3 %f -n 1 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> ref2.txt", thetanuc, SizeRatioC, SizeRatioA, coalT1, coalT2, coalT3))
  # locus 3
  system(sprintf("./ms 80 1 -t %f -I 4 40 30 8 2 -n 3 %f -n 1 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> ref2.txt", thetanuc, SizeRatioC, SizeRatioA, coalT1, coalT2, coalT3))
  # mitochondrial
  system(sprintf("./ms 60 1 -t %f -I 4 35 20 3 2 -n 3 %f -n 1 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2 | ./sample_stats >> ref2_mit.txt", thetamt, SizeRatioC, SizeRatioA, coalT1, coalT2, coalT3))
  
  ## save parameter values
  parameters <- rbind(parameters, data.frame(Ne, thetanuc, thetamt, coalT1, coalT2, coalT3))
}

## calculate summary statistics  
vic <- read.table("vic.txt")
vic_mit <- read.table("vic_mit.txt")
vic <- data.frame(pi.m=tapply(vic[,2]/Lnu, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  ss.m=tapply(vic[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  D.m=tapply(vic[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  T.m=tapply(vic[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  TH.m=tapply(vic[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  pi.v=tapply(vic[,2]/Lnu, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  ss.v=tapply(vic[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  D.v=tapply(vic[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  pi.m=vic_mit[,2]/Lmt,
                  ss.m=vic_mit[,4],
                  D.m=vic_mit[,6],
                  T.m=vic_mit[,8],
                  TH.m=vic_mit[,10])

p_d <- read.table("p_d.txt")
p_d_mit <- read.table("p_d_mit.txt")
p_d <- data.frame(pi.m=tapply(p_d[,2]/Lnu, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  ss.m=tapply(p_d[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  D.m=tapply(p_d[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  T.m=tapply(p_d[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  TH.m=tapply(p_d[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  pi.v=tapply(p_d[,2]/Lnu, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  ss.v=tapply(p_d[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  D.v=tapply(p_d[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                  pi.m=p_d_mit[,2]/Lmt,
                  ss.m=p_d_mit[,4],
                  D.m=p_d_mit[,6],
                  T.m=p_d_mit[,8],
                  TH.m=p_d_mit[,10])

ref1 <- read.table("ref1.txt")
ref1_mit <- read.table("ref1_mit.txt")
ref1 <- data.frame(pi.m=tapply(ref1[,2]/Lnu, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   ss.m=tapply(ref1[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   D.m=tapply(ref1[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   T.m=tapply(ref1[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   TH.m=tapply(ref1[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   pi.v=tapply(ref1[,2]/Lnu, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   ss.v=tapply(ref1[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   D.v=tapply(ref1[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   pi.m=ref1_mit[,2]/Lmt,
                   ss.m=ref1_mit[,4],
                   D.m=ref1_mit[,6],
                   T.m=ref1_mit[,8],
                   TH.m=ref1_mit[,10])

ref2 <- read.table("ref2.txt")
ref2_mit <- read.table("ref2_mit.txt")
ref2 <- data.frame(pi.m=tapply(ref2[,2]/Lnu, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   ss.m=tapply(ref2[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   D.m=tapply(ref2[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   T.m=tapply(ref2[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   TH.m=tapply(ref2[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   pi.v=tapply(ref2[,2]/Lnu, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   ss.v=tapply(ref2[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   D.v=tapply(ref2[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                   pi.m=ref2_mit[,2]/Lmt,
                   ss.m=ref2_mit[,4],
                   D.m=ref2_mit[,6],
                   T.m=ref2_mit[,8],
                   TH.m=ref2_mit[,10])

models <- rep(c("vic", "p_d", "ref1", "ref2"), each=numsim)

## join data
sust <- rbind(vic, p_d, ref1, ref2)
names(sust) <- c("pi.m", "ss.m", "D.m","T.m", "TH.m","pi.v","ss.v","D.v","pi.m.1", "ss.m.1", "D.m.1","T.m.1", "TH.m.1")
## load empirical data
# bfib
l1<-as.vector(c(5.124574,46,-1.232949,0.56648,4.558091))
names(l1) <- c("pi.m", "ss.m", "D.m", "T.m", "TH.m")
# 27_28
l2<-as.vector(c(4.850869,29,-0.343913,2.41039, 2.440476 ))
names(l2) <- c("pi.m", "ss.m", "D.m", "T.m", "TH.m")
# rpl3
l3<-as.vector(c(11.371835, 69,-0.609185, 8.24841, 3.123418 ))
names(l3) <- c("pi.m", "ss.m", "D.m", "T.m", "TH.m")
l<-rbind(l1,l2,l3)
# mitochondrial
mit<-as.vector(c(36.266102, 133, 0.951737, 57.327119, -21.061017 ))
names(mit) <- c("pi.m", "ss.m", "D.m", "T.m", "TH.m")

emp<- data.frame(pi.m=tapply(l[,1], factor(rep(1:1, each=numloc)), mean,na.rm=T),
                 ss.m=tapply(l[,2], factor(rep(1:1, each=numloc)), mean,na.rm=T),
                 D.m=tapply(l[,3], factor(rep(1:1, each=numloc)), mean,na.rm=T),
                 T.m=tapply(l[,4], factor(rep(1:1, each=numloc)), mean,na.rm=T),
                 TH.m=tapply(l[,5], factor(rep(1:1, each=numloc)), mean,na.rm=T),
                 pi.v=tapply(l[,1], factor(rep(1:1, each=numloc)), var,na.rm=T),
                 ss.v=tapply(l[,2], factor(rep(1:1, each=numloc)), var,na.rm=T),
                 D.v=tapply(l[,3], factor(rep(1:1, each=numloc)), var,na.rm=T),
                 pi.m=mit[1],
                 ss.m=mit[2],
                 D.m=mit[3],
                 T.m=mit[4],
                 TH.m=mit[5])
names(emp) <- c("pi.m", "ss.m", "D.m","T.m", "TH.m","pi.v","ss.v","D.v","pi.m.1", "ss.m.1", "D.m.1","T.m.1", "TH.m.1")

# save models, parameters and sust
write.table(models, file="models.txt", quote=F, row.names=F, col.names=F)
write.table(parameters, file="parameters.txt", quote=F, row.names=F, col.names=F)
write.table(sust, file="sust.txt", quote=F, row.names=F, col.names=F)

# cross-validation test with untransformed and transformed by Principal Component Analysis (PCA) datasets to assess their performance 
# transformation of summary statistics (sust) by PCA 
pcasust <- prcomp(sust, scale = TRUE)
pcaemp<-predict(pcasust, emp, scale=TRUE)
summary(pcasust) # check each PC importance
pcasust<-pcasust$x[,1:6] # select the number of PCs 
pcaemp<-pcaemp[,1:6] # repeat the number of PCs for the empirical data

# load abc package
library(abc)
 
# untransformed dataset
cvlog <- cv4postpr(models, sust, nval=10, tol= 0.05, method="mnlogistic")
summary(cvlog)

# transformed dataset
cvlogpca <- cv4postpr(models, pcasust, nval=10, tol 0.05), method="mnlogistic")
summary(cvlogpca)

# model selection using the algorithm with best performance according to cross-validation test 

LOG<-postpr(emp, models, sust, tol = 0.05, method = "mnlogistic")
PCALOG<-postpr(pcaemp, models, pcasust, tol = 0.05, method = "mnlogistic")




