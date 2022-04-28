#!/bin/env Rscript
# author: ph-u
# script: interactionTypes.r
# desc: interaction types through time
# in: Rscript interactionTypes.r [path/to/data/] [data_basename] [math-model-type] [num_accepted_simulations]
# out: result/*-eco*.csv
# arg: 3
# date: 20220109

argv=(commandArgs(T))
source("../pipeline/src.r")
library(deSolve)
pT = argv[1]; nAm = argv[2]; tYpe=argv[3]
p0 = read.csv(paste0(pT,nAm,"-sam.csv"), header=T, stringsAsFactors=F)
pR = read.csv(paste0(pT,nAm,"-pri.csv"), header=T, stringsAsFactors=F)
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T)
nAcc = min(as.numeric(argv[4]), nrow(p0))

##### trim illogical "best parameters" ##### 20220411
if(tYpe=="LVC"){oDe="c"}else{oDe="g"}
x0 = rep(0,ncol(t0)-1)
for(i in 2:ncol(t0)){x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])}
pRM = c();for(i in 1:nrow(p0)){
	a0 = solveLV(x0, as.numeric(p0[i,]), range(t0[,1]), oDe)
	for(i0 in 2:ncol(a0)){a0[,i0] = ifelse(a0[,i0]>150 | a0[,i0]<0,-100,a0[,i0])};rm(i0)
	a0[is.na(a0)] = -100
	if(any(a0==-100)){pRM = c(pRM,i)}
};rm(i)
if(length(pRM)>0){p0=p0[-pRM,]}
p0 = p0[,which(pR[,2] != "r" & pR[,2] != "k")]
pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]

##### single ecological relationship ##### 2022{0109,0111,0331}
tY = c("competition","mutualism","neutral/no_interaction","predator/parasite_of_other","prey/host_of_other","commensal_of_other","commensal_Host_of_other","harmed_by_other","harming_other")
sEr = function(i1,i2,cT=tY,tP=tYpe){
	if(tP=="LVC"){ # LVC: + harm ; - help ; compare i1 (1->2) to i2 (2->1)
		iT = ifelse(i1<0, ifelse(i2<0,tY[2], ifelse(i2>0,tY[5],tY[7])),
		ifelse(i1>0, ifelse(i2>0,tY[1],ifelse(i2<0,tY[4],tY[9])),
		ifelse(i2>0,tY[8],ifelse(i2<0,tY[6],tY[3]))))
	}else{ # gLV: - harm ; + help ; compare i1 (1->2) to i2 (2->1)
		iT = ifelse(i1>0, ifelse(i2>0,tY[2], ifelse(i2<0,tY[5],tY[7])),
		ifelse(i1<0, ifelse(i2<0,tY[1],ifelse(i2>0,tY[4],tY[9])),
		ifelse(i2<0,tY[8],ifelse(i2>0,tY[6],tY[3]))))
	}
	return(iT)
}

##### map columns into interaction matrix #####
n = unique(pR[,1])
if(tYpe=="LVC"){
	a = as.data.frame(matrix(1, nr=length(n), nc=length(n)))
	i0=1
	for(i in 1:nrow(a)){ for(j in 1:ncol(a)){if(i!=j){
		a[i,j] = colnames(p0)[i0]
		i0=i0+1
}}}}else{
	a = as.data.frame(t(matrix(colnames(p0), nr=length(n), nc=length(n))))
};row.names(a) = colnames(a) = n

##### proportion of deduced ecological interaction types #####
w = as.data.frame(matrix(0, nr=ifelse(tYpe=="LVC",ncol(p0)/2,length(n)*(length(n)-1)), nc=length(tY)+4))
colnames(w) = c("subject","other",tY,"Chi-sq","p-val")
if(tYpe=="LVC"){
	i9 = 2:length(n)
	j9 = 1:(length(n)-1)
}else{i9 = j9 = 1:length(n)}
i0=1; for(i in i9){for(j in j9){if(ifelse(tYpe=="LVC",i>j,i>=j)){
	w[i0,1:2] = c(n[i],n[j])
	b = table(sEr(p0[,a[i,j]],p0[,a[j,i]]))
	w[i0,names(b)] = b
	if(length(b)>1){
		l = chisq.test(b)
		w[i0,(ncol(w)-1):ncol(w)] = c(l$statistic,l$p.value)
		if(length(b)>2 && l$p.value<.05){ ## pairwise Chi-sq tests with bonferroni p-value correction
			e = as.data.frame(matrix(NA, nr=choose(length(b),2), nc=4))
			colnames(e) = c("interaction1","interaction2",rev(rev(colnames(w))[1:2]))
			j0=1; for(i1 in 1:(length(b)-1)){for(j1 in 2:length(b)){if(i1<j1){
				b0 = chisq.test(b[c(i1,j1)])
				e[j0,] = c(names(b)[c(i1,j1)],b0$statistic,b0$p.value*nrow(e))
				j0=j0+1}}}
			write.csv(e, paste0(pT,"../result/",nAm,"-eco_",n[i],"_",n[j],".csv"), quote=F, row.names=F)
	}}else{w[i0,(ncol(w)-1):ncol(w)] = rep(NA,2)};i0=i0+1}}}
write.csv(w, paste0(pT,"../result/",nAm,"-eco.csv"), quote=F, row.names=F)
