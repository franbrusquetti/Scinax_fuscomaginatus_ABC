## in order to use this code you need to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

## set working directory  
setwd(“home/Documents/Simulation_Data”)

# Goodness-of-fit (gfit) of the data simulated from the best-supported model (model 3: isolation by instability 1) in relation to empirical data
# load summary statistics (sust) and models 
sust <- read.table("sust.txt")
models <- read.table("models.txt")

# extract summary statistics of ref1 model
sust_ref1 <- subset(sust, subset=models=="ref1")

# plot to check the adequacy of the results with the values of the observed data within the values of the simulated data
# gfit with untransformed dataset
gfit_ref1<-gfit(target=emp, sumstat=sust_ref1, nb.replicate=1000, tol=0.05, statistic=mean)
summary(gfit_ref1)
plot(gfit_ref1)

# gfit with transformed by PCA dataset
pcasust_ref1<- prcomp(sust_ref1, scale = TRUE)
pcaemp_ref1<-predict(pcasust_ref1, emp, scale=TRUE) 
summary(pcasust_ref1) # check each PC importance
pcasust_ref1<-pcasust_ref1$x[,1:6] # select the number of PCs 

pcaemp_ref1<-pcaemp_ref1[,1:6] # repeat the number of PCs for the empirical data

gfit_pcaref1<-gfit(target=pcaemp_ref1, sumstat=pcasust_ref1, nb.replicate=100, tol=0.05, statistic=mean)
summary(gfit_pcaref1)
plot(gfit_pcaref1)

# Posterior predictive check (PPC) of each summary statistic simulated from the best-supported model (model 3: isolation by instability 1) in relation to empirical data

### variable declarations

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

## number of simulations for the posterior check
simpost<-1000

## load parameters of the selected model and empirical data
parameters<-read.table("parameters.txt")
parameters_ref1<-subset(parameters, subset=models=="ref1")
names(parameters_ref1)<- c("Ne", "thetanuc", "thetamt", "SizeRatioC", "SizeRatioA", "SizeRatioD", "SizeRatioB","coalT1", "coalT2", "coalT3")
names(sust_ref1) <- c("pi.m", "ss.m", "D.m","T.m", "TH.m","pi.v","ss.v","D.v","pi.m.1", "ss.m.1", "D.m.1","T.m.1", "TH.m.1")

## load empirical data (emp; see the script 1_runms_final)

## perform abc for inferring the parameters
peNN<-abc(target=emp, param=parameters_ref1, sumstat=sust_ref1, tol=0.05, method="neuralnet",transf=c("logit","logit","logit","logit", "logit","logit", "logit","logit", "logit", "logit"),logit.bounds = rbind(c(2477996, 5976161),c(19.15,64.95),c(31.74,121.15),c(0.16,0.95), c(0.16, 0.95),c(0.16,0.95), c(0.16, 0.95), c(0.002, 0.72), c(0.009, 1.05), c(0.02, 1.51)))

## sampling with replacement in the multivariate posterior distribution
newsamp<-sample(1:(dim(peNN$adj)[1]),size=simpost,replace=T,prob=peNN$weights)

newsamp<-peNN$adj[newsamp,]

ppc <- data.frame()

## Counter
linecounter = 1

for (i in 1:simpost) {
 
  ### Define parameters
  thetanuc <- newsamp[linecounter,2]
  thetamt <-newsamp[linecounter,3]
  SizeRatioC<-newsamp[linecounter,4]
  SizeRatioA<-newsamp[linecounter,5]
  SizeRatioD<-newsamp[linecounter,6]
  SizeRatioB<-newsamp[linecounter,7]
  coalT1 <- newsamp[linecounter,8]/newsamp[linecounter,1]
  coalT2 <- newsamp[linecounter,9]/newsamp[linecounter,1]
  coalT3<-newsamp[linecounter,10]/newsamp[linecounter,1]

    
  ## ms's command: 4 populations with bottleneck in all populations, no migration
  # locus 1
  system(sprintf("./ms 124 1 -t %f -I 4 72 42 6 4 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ppc_nuc.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # locus 2
  system(sprintf("./ms 114 1 -t %f -I 4 72 36 2 4 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ppc_nuc.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # locus 3
  system(sprintf("./ms 80 1 -t %f -I 4 40 30 8 2 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ppc_nuc.txt", thetanuc, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
  # mitochondrial
  system(sprintf("./ms 60 1 -t %f -I 4 35 20 3 2 -n 3 %f -n 1 %f -n 4 %f -n 2 %f -ej %f 3 1 -ej %f 1 4 -ej %f 4 2  | ./sample_stats >> ppc_mit.txt", thetamt, SizeRatioC, SizeRatioA, SizeRatioD, SizeRatioB, coalT1, coalT2, coalT3))
 
	linecounter <- linecounter +1
	}

# calculate the summary statistics
ppc_nuc <- read.table("ppc_nuc.txt")
ppc_mit <- read.table("ppc_mit.txt")
ppc <- data.frame(pi.m=tapply(ppc_nuc[,2]/Lnu, factor(rep(1:simpost, each=numloc)), mean,na.rm=T),
              	ss.m=tapply(ppc_nuc[,4], factor(rep(1:simpost, each=numloc)), mean,na.rm=T),
              	D.m=tapply(ppc_nuc[,6], factor(rep(1:simpost, each=numloc)), mean,na.rm=T),
              	T.m=tapply(ppc_nuc[,8], factor(rep(1:simpost, each=numloc)), mean,na.rm=T),
              	TH.m=tapply(ppc_nuc[,10], factor(rep(1:simpost, each=numloc)), mean,na.rm=T),
              	pi.v=tapply(ppc_nuc[,2]/Lnu, factor(rep(1:simpost, each=numloc)), var,na.rm=T),
              	ss.v=tapply(ppc_nuc[,4], factor(rep(1:simpost, each=numloc)), var,na.rm=T),
              	D.v=tapply(ppc_nuc[,6], factor(rep(1:simpost, each=numloc)), var,na.rm=T),
              	pi.m.1=ppc_mit[,2]/Lmt,
              	ss.m.1=ppc_mit[,4],
              	D.m.1=ppc_mit[,6],
              	T.m.1=ppc_mit[,8],
              	TH.m.1=ppc_mit[,10])

# plot the previously selected variables

hist(ppc[,2],breaks=40, xlab="ss.m", main=""); abline(v = emp[,1], col = 2)
hist(ppc[,4],breaks=40, xlab=" T.m ", main=""); abline(v = emp[,2], col = 2)
hist(ppc[,5],breaks=40, xlab=" TH.m", main=""); abline(v = emp[,3], col = 2)
hist(ppc[,7],breaks=40, xlab=" ss.v", main=""); abline(v = emp[,4], col = 2)
hist(ppc[,10],breaks=40, xlab="ss.m.1  ", main=""); abline(v = emp[,5], col = 2)
hist(ppc[,12],breaks=40, xlab="T.m.1", main=""); abline(v = emp[,6], col = 2)
hist(ppc[,13],breaks=40, xlab=" TH.m.1", main=""); abline(v = emp[,7], col = 2)

