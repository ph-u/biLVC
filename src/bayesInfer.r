#!/bin/env Rscript
# author: ph-u
# script: bayesInfer.r
# desc: Bayesian Inference on Lotka-Volterra Competition/generalized Model
# in: Rscript bayesInfer.r [basename] [competition (c) / generalized (g)] [max-iteration (in 10^5)] [replicate No.] [(clusster only) /full/path]
# out: result/*.pdf, data/*.{csv,txt}
# arg: 4/5
# date: 20220107

argv=(commandArgs(T))
kk0 = read.csv("../parIN.csv", header = F)
kk1 = as.numeric(argv[1])
library(FME) # FME(1.3.6.2), deSolve(1.30), rootSolve(1.8.2.3), coda(0.19.4)
# doi: 10.18637/jss.v033.i03
pT=paste0(ifelse(version$os=="linux-gnu",paste0(kk0[kk1,5],"/"),""),"../")
source(paste0(pT,"src/src.r"))
nM=kk0[kk1,1]; tP=kk0[kk1,2]
sT=read.csv(paste0(pT,"data/",nM,"-log.csv"), header=T)
rEc=read.csv(paste0(pT,"data/",nM,"-pri.csv"), header=T)
sEed=read.csv(paste0(pT,"data/",nM,"-seed.csv"), header=T)
set.seed(sEed$seed[kk0[kk1,4]])

##### MCMC ##### 20220108, 20220331 (+oDe option)
mcMC = modMCMC(f=mCres, rEc[,"initial"], df=sT, oDe=tP, lower=rEc[,"min"], upper=rEc[,"max"], niter=kk0[kk1,3]*10^5, outputlength=max(1e2,kk0[kk1,3]*10^2), updatecov=50, burninlength=0)

##### summary ##### 20220108
#save(mcMC, file=paste0(pT,"data/",nM,"-",kk0[kk1,4],"-sam.RData"))
write.csv(as.data.frame(mcMC$par), paste0(pT,"data/",nM,"-",kk0[kk1,4],"-sam.csv"), quote=F, row.names=F)
write.csv(summary(mcMC), paste0(pT,"result/",nM,"-",kk0[kk1,4],"-est.csv"), quote=F, row.names=T)

pdf(paste0(pT,"result/",nM,"-",kk0[kk1,4],"-sum.pdf"))
plot(mcMC)
invisible(dev.off())

pdf(paste0(pT,"result/",nM,"-",kk0[kk1,4],"-hist.pdf"))
hist(mcMC, breaks=20)
invisible(dev.off())

pdf(paste0(pT,"result/",nM,"-",kk0[kk1,4],"-cumu.pdf"))
cumuplot(as.mcmc(mcMC$pars))
invisible(dev.off())

