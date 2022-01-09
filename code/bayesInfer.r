#!/bin/env Rscript
# author: ph-u
# script: bayesInfer.r
# desc: Bayesian Inference on Lotka-Volterra Competition Model
# in: Rscript bayesInfer.r [basename]
# out: result/*.pdf, data/*.{csv,txt}
# arg: 1
# date: 20220107

#SBATCH -J mcmcR-C
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH -p skylake

argv=(commandArgs(T))
library(FME) # FME(1.3.6.2), deSolve(1.30), rootSolve(1.8.2.3), coda(0.19.4)
# doi: 10.18637/jss.v033.i03
pT=ifelse(version$os=="linux-gnu","/home/pmh65/rds/0_abcrPar/","../")
source(paste0(pT,"code/src.r"))
nM=argv[1]
sT=read.csv(paste0(pT,"raw/",nM,".csv"), header=T)

##### log-transformation ##### 20220101
sT[,-1] = log(sT[,-1])
write.csv(sT,paste0(pT,"data/",nM,"-log.csv"), quote=F, row.names=F)

##### Rough priors ##### 20220101
nSp = ncol(sT)-1
rEc = as.data.frame(matrix(NA,nr=nSp*(nSp+2),nc=5))
colnames(rEc) = c("id","influencer","initial","min","max")
rEc[,1] = rep(colnames(sT)[-1], each=nSp+2)
rEc[,2] = rep(c("r","k",colnames(sT)[-1]),nSp)
rEc = rEc[which(rEc[,1]!=rEc[,2]),] ## fix self-interaction =1
for(i in 1:nrow(rEc)){
        if(rEc[i,2]=="r"){
                rEc[i,3:ncol(rEc)] = rolRate(sT,rEc[i,1],1)
        }else if(rEc[i,2]=="k"){
                p = range(sT[,rEc[i,1]])
                rEc[i,3:ncol(rEc)] = c(mean(p),0,p[2]+.2*mean(p))
        }else{
                rEc[i,3:ncol(rEc)] = c(0,-2,2)
        }
}

##### Fine parameter estimation ##### 20220108
mcFIT = modFit(oBj,rEc[,3], sT, bestMet(rEc[,3], sT), method="Nelder-Mead")
rEc[,"initial"] = mcFIT$par
rEc[,"min"] = ifelse(rEc[,"min"]>=mcFIT$par, mcFIT$par-abs(rEc[,"min"]), rEc[,"min"])
rEc[,"max"] = ifelse(rEc[,"max"]<=mcFIT$par, mcFIT$par+abs(rEc[,"max"]), rEc[,"max"])
write.csv(rEc, paste0(pT,"data/",argv[1],"-pri.csv"), quote=F, row.names=F)

##### MCMC ##### 20220108
mcMC = modMCMC(f=mCres, rEc[,"initial"], df=sT, lower=rEc[,"min"], upper=rEc[,"max"], niter=1e5, outputlength=1e2, updatecov=50, burninlength=0)

##### summary ##### 20220108
pdf(paste0(pT,"result/",argv[1],"-sum.pdf"))
plot(mcMC)
invisible(dev.off())

pdf(paste0(pT,"result/",argv[1],"-hist.pdf"))
hist(mcMC, breaks=20)
invisible(dev.off())

write.table(summary(mcMC), paste0(pT,"data/",argv[1],"-mcmcSum.txt"), quote=F, row.names=F, sep=",")

pdf(paste0(pT,"result/",argv[1],"-cumu.pdf"))
cumplot(as.mcmc(mcMC$pars))
invisible(dev.off())

