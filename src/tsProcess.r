#!/bin/env Rscript
# author: ph-u
# script: tsProcess.r
# desc: time-series raw data process
# in: Rscript tsProcess.r [basename] [competition (c) / generalized (g)] [replicates] [/full/path (optional)] [raw ts directory name]
# out: data/[basename]-{log,pri,seed}.csv
# arg: 3/4
# date: 20220508 (segregate from bayesInfer.r)

argv=(commandArgs(T))
kk0 = read.csv("../parIN.csv", header = F)
kk1 = as.numeric(argv[1])
library(FME)
pT=paste0(ifelse(version$os=="linux-gnu",paste0(kk0[kk1,4],"/"),""),"../")
source(paste0(pT,"src/src.r"))
nM=kk0[kk1,1]; tP=kk0[kk1,2]
sT=read.csv(paste0(pT,kk0[kk1,5],"/",nM,".csv"), header=T)

##### transformation ##### 20220101, 20220331 (+ % data type)
if(all(sT[,-1]<=1)){sT[,-1] = sT[,-1] * 100}else{sT[,-1] = log(sT[,-1]+1)}
write.csv(sT,paste0(pT,"data/",nM,"-log.csv"), quote=F, row.names=F)

##### Rough priors ##### 20220101, 20220331 (+gLV parameter version)
nSp = ncol(sT)-1
rEc = as.data.frame(matrix(NA,nr=nSp*(nSp+ifelse(tP=="c",2,1)),nc=5))
colnames(rEc) = c("id","influencer","initial","min","max")
rEc[,1] = rep(colnames(sT)[-1], each=nSp+ifelse(tP=="c",2,1))
if(tP=="c"){
        rEc[,2] = rep(c("r","k",colnames(sT)[-1]),nSp)
        rEc = rEc[which(rEc[,1]!=rEc[,2]),] ## fix self-interaction =1
}else{
        rEc[,2] = rep(c("r",colnames(sT)[-1]),nSp)
}
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

##### Fine parameter estimation ##### 20220108, 20220331 (+oDe option)
mcFIT = modFit(oBj,rEc[,3], sT, tP, bestMet(rEc[,3], sT, tP), method="Nelder-Mead")
rEc[,"initial"] = mcFIT$par
rEc[,"min"] = ifelse(rEc[,"min"]>=mcFIT$par, mcFIT$par-abs(rEc[,"min"]), rEc[,"min"])
rEc[,"max"] = ifelse(rEc[,"max"]<=mcFIT$par, mcFIT$par+abs(rEc[,"max"]), rEc[,"max"])
write.csv(rEc, paste0(pT,"data/",nM,"-pri.csv"), quote=F, row.names=F)

##### set seeds ##### 20220509
sEed = data.frame("replicate"=1:kk0[kk1,3],"seed"=ceiling(runif(kk0[kk1,3])*10^9))
write.csv(sEed,paste0(pT,"data/",nM,"-seed.csv"), quote=F, row.names=F)
