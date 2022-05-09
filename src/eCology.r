#!/bin/env Rscript
# author: ph-u
# script: eCology.r
# desc: extract ecology from simulation replicates
# in: Rscript eCology.r [path/to/data/] [time_series_basename] [LVC/gLV] [replicates]
# out: result/*-eco.csv
# arg: 4
# date: 20220508 (supersede interactionTypes.r)

##### env #####
argv=(commandArgs(T))
source("src.r")
library(deSolve)
library(PMCMRplus) # v1.9.3
pT = argv[1]; nAm = argv[2]; tYpe=argv[3]
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T)
pR = read.csv(paste0(pT,nAm,"-pri.csv"), header=T)
rP = vector(mode="list", length=as.numeric(argv[4]))
for(i in 1:as.numeric(argv[4])){
	rP[[i]] = read.csv(paste0(pT,nAm,"-",i,"-sam.csv"), header=T)
}

##### categorize ecology ##### (from interactionTypes.r)
tY = c("mutualism","commensal_Host_of_c2","prey/host_of_c2","commensal_of_c2","neutral/no_interaction","harmed_by_c2","predator/parasite_of_c2","harming_c2","competition")
sEr = function(i1,i2,cT=tY,tP=tYpe){
        rEf = data.frame("P2p"=rep(1:-1,each=3),"p2P"=rep(1:-1,3),"P_is_the"=cT)
        if(tP=="LVC"){i1 = -i1;i2 = -i2}
        i1 = ifelse(i1==0,0,ifelse(i1>0,1,-1))
        i2 = ifelse(i2==0,0,ifelse(i2>0,1,-1))
        iT = rep(NA,length(i1))
        for(i in 1:length(iT)){iT[i] = rEf$P_is_the[which(rEf$P2p==i1[i] & rEf$p2P==i2[i])]}
        return(iT)
}

##### exclude non-fitting simulations ##### 20220411
oDe = ifelse(tYpe=="LVC","c","g")
x0 = rep(0,ncol(t0)-1)
for(i in 2:ncol(t0)){x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])}
for(i1 in 1:length(rP)){
	pRM = c();for(i in 1:nrow(rP[[i1]])){
		a0 = solveLV(x0, as.numeric(rP[[i1]][i,]), range(t0[,1]), oDe)
	        for(i0 in 2:ncol(a0)){a0[,i0] = ifelse(a0[,i0]>150 | a0[,i0]<0,-100,a0[,i0])};rm(i0)
		a0[is.na(a0)] = -100
		if(any(a0==-100)){pRM = c(pRM,i)}
	};if(length(pRM)>0){rP[[i1]]=rP[[i1]][-pRM,]}
	rP[[i1]] = rP[[i1]][,which(pR[,2] != "r" & pR[,2] != "k")]
};rm(i,i1)
pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]

##### interaction matrix ##### (from interactionTypes.r)
n = unique(pR[,1])
if(tYpe=="LVC"){
        a = as.data.frame(matrix(1, nr=length(n), nc=length(n)))
        i0=1
        for(i in 1:nrow(a)){ for(j in 1:ncol(a)){if(i!=j){
                a[i,j] = colnames(rP[[1]])[i0]
                i0=i0+1
}}}}else{
        a = as.data.frame(t(matrix(colnames(rP[[1]]), nr=length(n), nc=length(n))))
};row.names(a) = colnames(a) = n

##### category pairwise combinations #####
catComb = as.data.frame(matrix(NA,nr=choose(length(n),2)+ifelse(tYpe=="LVC",0,length(n)), nc=2))
nR=i=j=1;repeat{
	if(ifelse(tYpe=="LVC",i<j,i<=j)){catComb[nR,] = c(n[i],n[j]);nR = nR+1}
	j = j+1
	if(j>length(n)){j=1;i = i+1}
	if(nR>nrow(catComb)){break}
};rm(nR,i,j)

##### interaction categorization result collector #####
eCo = as.data.frame(matrix(0,nr=(choose(length(n),2)+ifelse(tYpe=="LVC",0,length(n)))*length(rP)*length(tY),nc=7))
colnames(eCo) = c(paste0("category",1:2),"replicate","c1_is","count","fit_sim","ratio_in_rep")
eCo$replicate = rep(1:length(rP),each=nrow(eCo)/length(rP))
eCo$c1_is = rep(tY,nrow(eCo)/length(tY))
eCo$category1 = rep(rep(catComb[,1],each=length(tY)),length(rP))
eCo$category2 = rep(rep(catComb[,2],each=length(tY)),length(rP))

##### map interactions #####
for(i in 1:length(rP)){
	eCo$fit_sim[which(eCo$replicate==i)] = nrow(rP[[i]])
	for(c2 in 1:length(n)){ for(c1 in 1:length(n)){
		if(ifelse(tYpe=="LVC",c2>c1,c2>=c1)){
			cAt = table(sEr(rP[[i]][,a[c2,c1]],rP[[i]][,a[c1,c2]])) # a[row,col] - "col" (P) affect population of "row" (p)
			for(i0 in 1:length(cAt)){
				eCo$count[which(eCo$category1==colnames(a)[c1] & eCo$category2==colnames(a)[c2] & eCo$replicate==i & eCo$c1_is==names(cAt)[i0])] = cAt[i0]
}
		}
	}}
};rm(i)
eCo$ratio_in_rep = eCo$count/eCo$fit_sim # relationship ratio in the top 100 best-fit after double biological simulation-data reality check
write.csv(eCo,paste0(pT,"../result/",nAm,"-eco.csv"), quote=F, row.names=F)

##### Kruskal test + posthoc Nemenyi (single-step p-adj) #####
for(i in 1:nrow(catComb)){
	i0 = eCo[which(eCo$category1==catComb[i,1] & eCo$category2==catComb[i,2]),]
	kW = kwAllPairsNemenyiTest(i0$ratio_in_rep~as.factor(i0$c1_is), dist="Chisquare")
	kwSta = as.data.frame(kW$statistic)
	kwPva = as.data.frame(kW$p.value)
	kwS = kwP = c();for(j in 1:ncol(kwSta)){
		kwS = c(kwS,kwSta[,j])
		kwP = c(kwP,kwPva[,j])
	};kwTab = data.frame(
		"interaction1"=rep(row.names(kwSta),ncol(kwSta)),
		"interaction2"=rep(colnames(kwSta),each=nrow(kwSta)),
		"chi.sq"=kwS, "adj.p"=kwP)
	kwTab = kwTab[which(!is.na(kwTab$adj.p) & kwTab$adj.p<=.1),]
	write.csv(kwTab,paste0(pT,"../result/",nAm,"-kwPairs_",catComb[i,1],"_",catComb[i,2],".csv"), quote=F, row.names=F)

## grouped boxplot
	pdf(paste0(pT,"../result/",nAm,"-kwPairs_",catComb[i,1],"_",catComb[i,2],".pdf"))
	par(mar=c(10,4.5,0,2)+.1, xpd=T)
	boxplot(i0$ratio_in_rep~gsub("/"," / ",gsub("_"," ",gsub("_c2","",gsub("_of_c2","",i0$c1_is)))),
	ylim=c(0,1+nrow(kwTab)/10), col="#FFFFFFFF", xlab="", ylab=paste0("ratio of likeliness in ",argv[4]," replicates"), las=2, pch=4, yaxt="n")
	axis(2,at=seq(0,1,.2),labels=seq(0,1,.2))
	lB = tY[order(tY)]
	lBplt = kwTab
	for(i2 in 1:length(lB)){for(i1 in 1:2){lBplt[which(lBplt[,i1]==lB[i2]),i1]=i2}}
	lBplt$chi.sq = round(lBplt$chi.sq,2)
	lBplt$adj.p = ifelse(lBplt$adj.p<.001,"<<0.01",round(lBplt$adj.p,3))
	segments(x0=as.numeric(lBplt[,1]),x1=as.numeric(lBplt[,2]),y0=(1:nrow(lBplt))/10+1)
	text(x = (as.numeric(lBplt[,1])+as.numeric(lBplt[,2]))/2-.5, y = (1:nrow(lBplt))/10+1.03, labels=paste("X =",lBplt[,3],"; adj-p =",lBplt[,4]), xpd=T)
	invisible(dev.off())
}
