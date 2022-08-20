#!/bin/env Rscript
# author: ph-u
# script: eCology.r
# desc: extract ecology from simulation replicates
# in: Rscript eCology.r [path/to/data/] [time_series_basename] [LVC/gLV] [replicates]
# out: result/*-{eco,tsMatch,kwPairs_[c1]_[c2]}.csv, result/*-{kwPairs_[c1]_[c2],tsAllRep}.pdf
# arg: 4
# date: 20220508 (supersede interactionTypes.r), 20220820 (CSD3 adaptation)

#SBATCH -J ecolAna
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p skylake-himem

##### env #####
argv=(commandArgs(T))
source("src.r")
library(deSolve)
#library(PMCMRplus) # v1.9.3
#source("fdrBH.r")
pT = argv[1]; nAm = argv[2]; tYpe=argv[3]; pOut = gsub("data","result",pT)
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T)
pR = read.csv(paste0(pT,nAm,"-pri.csv"), header=T)
sD = read.csv(paste0(pT,nAm,"-seed.csv"), header=T)
rP = rK = vector(mode="list", length=as.numeric(argv[4]))
for(i in 1:as.numeric(argv[4])){
	rP[[i]] = read.csv(paste0(pT,nAm,"-",i,"-sam.csv"), header=T)
}
dMaxSim = nrow(rP[[1]])*as.numeric(argv[4]) # ref max simulation

##### colour ##### (from SimDataPlot.r)
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)
cBl = palette.colors(palette = "Okabe-Ito", alpha=.1, recycle = T)
n=colnames(t0)[-1]
cOl = data.frame(id=n,cPt=cBp[1:length(n)],cLn=cBl[1:length(n)])
ptCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),2]
lnCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),3]
if(any(t0[,-1]>30)){yLab="percentage presence [%]"}else{yLab="log_e(y+1) [CFU/mL]"}
if(argv[4]=="LVC"){oDe="c"}else{oDe="g"}
ptCol0 = rep(ptCol,ceiling((ncol(t0)-1)/length(ptCol)))
lnCol0 = rep(lnCol,ceiling((ncol(t0)-1)/length(lnCol)))

##### categorize ecology ##### (from interactionTypes.r)
tY0 = c("mutu","cHos","prey","comm","neut","hmED","pred","hmIG","comp")
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

##### get simulation match range ##### (from SimDataPlot.r)
acRatio = .1
tUq = unique(t0[,1])
tUq = tUq[order(tUq)[-1]] ## accending order safety net
dMin = dMax = dRec = as.data.frame(matrix(0,nr=length(tUq),nc=ncol(t0)))
colnames(dMin) = colnames(dMax) = colnames(dRec) = colnames(t0)
dMin[,1] = dMax[,1] = dRec[,1] = tUq
for(i in 1:length(tUq)){
        d = t0[which(t0[,1]==tUq[i]),]
        for(j in 2:ncol(t0)){
                d0 = range(d[,j])
                dMin[i,j] = ifelse(length(t0[,1])==length(unique(t0[,1])),min(d0[1]*(1-acRatio),d0[1]-1),max(d0[1]-diff(d0)/2,0))
                dMax[i,j] = ifelse(length(t0[,1])==length(unique(t0[,1])),max(d0[2]*(1+acRatio),d0[2]+1),d0[2]+diff(d0)/2)
}}

##### exclude non-fitting simulations ##### 20220411
oDe = ifelse(tYpe=="LVC","c","g")
x0 = rep(0,ncol(t0)-1)
for(i in 2:ncol(t0)){x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])}

##### plot Time-series #####
pdf(paste0(pOut,nAm,"-tsAllRep.pdf"), width=14)
par(mar=c(5,5,1,12)+.1, xpd=T)
matplot(t0[,1],t0[,-1], type="p", pch=(1:(ncol(t0)-1))%%25, cex=1.2, col=ptCol0,
        xlab=paste0(gsub("_"," (",colnames(t0)[1]),ifelse(length(grep("_",colnames(t0)[1]))>0,")","")),
        ylab=yLab, cex.axis=2, cex.lab=2)
legend("topright", inset=c(-.19,0), legend = colnames(t0)[-1], pch = (1:(ncol(t0)-1))%%25, lty=(1:(ncol(t0)-1))%%5+1, lwd=2, col = ptCol0)

i9=0;for(i1 in 1:length(rP)){
	set.seed(sD$seed[i1])
	pRM = c();for(i in 1:nrow(rP[[i1]])){
		a0 = solveLV(x0, as.numeric(rP[[i1]][i,]), range(t0[,1]), oDe)
	        for(i0 in 2:ncol(a0)){a0[,i0] = ifelse(a0[,i0]>150 | a0[,i0]<0,-100,a0[,i0])};rm(i0)
		a0[is.na(a0)] = -100
		if(any(a0==-100)){pRM = c(pRM,i)}else{
			matplot(a0[,1],a0[,-1], type="l", add=T, col=lnCol0, lty=(1:(ncol(t0)-1))%%5+1)
			i9=i9+1
		}
## Simulation-data match count
		for(j in 1:length(tUq)){for(k in 2:ncol(t0)){
			dRec[j,k] = dRec[j,k] + (a0[j,k]>=dMin[j,k] & a0[j,k]<=dMax[j,k])

        }}};if(length(pRM)>0){rP[[i1]]=rP[[i1]][-pRM,]}
	rK[[i1]] = rP[[i1]][,which(pR[,2] == "r" | pR[,2] == "k")]
	rP[[i1]] = rP[[i1]][,which(pR[,2] != "r" & pR[,2] != "k")]
};rm(i,i1)
pRrk = pR[which(pR[,2] == "r" | pR[,2] == "k"),]
pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]

##### ex: growth rate, carrying capacity #####
if(nrow(rK[[1]])>0){
	rkVal = rK[[1]]
	rkVal$replicate = 1
}else{
	rkVal = as.data.frame(matrix(0,nr=0,nc=nrow(pR)+1))
	colnames(rkVal)[ncol(rkVal)]="replicate"
}
for(i in 2:length(rK)){if(nrow(rK[[i]])>0){
	tMp = rK[[i]];tMp$replicate = i
	rkVal = rbind(rkVal,tMp)
}};colnames(rkVal)[-ncol(rkVal)] = paste(pRrk[,1],pRrk[,2],sep=".")
write.csv(rkVal,paste0(pT,nAm,"-rkValue.csv"), quote=F, row.names=F)

##### ex: simulation percentage match on data ##### (from SimDataPlot.r)
dRec[,-1] = dRec[,-1]/dMaxSim
write.csv(dRec,paste0(pT,nAm,"-tsMatch.csv"), quote=F, row.names=F)

##### ex: Time-series plot #####
text(min(t0[,1])+diff(range(t0[,1]))*.25,max(t0[,-1]),paste("Number of simulation(s)\nPlotted:",i9, "set(s)"), cex=1.2)
invisible(dev.off())

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
		if(ifelse(tYpe=="LVC",c2>c1,c2>=c1) & nrow(rP[[i]])>0){
			cAt = table(sEr(rP[[i]][,a[c2,c1]],rP[[i]][,a[c1,c2]])) # a[row,col] - "col" (P) affect population of "row" (p)
			for(i0 in 1:length(cAt)){
				eCo$count[which(eCo$category1==colnames(a)[c1] & eCo$category2==colnames(a)[c2] & eCo$replicate==i & eCo$c1_is==names(cAt)[i0])] = cAt[i0]
}
		}
	}}
};rm(i)
eCo$ratio_in_rep = ifelse(eCo$fit_sim==0,0,eCo$count/eCo$fit_sim) # relationship ratio in the top 100 best-fit after double biological simulation-data reality check
write.csv(eCo,paste0(pT,nAm,"-eco.csv"), quote=F, row.names=F)

##### Kruskal test + posthoc Nemenyi (single-step p-adj) #####
#for(i in 1:nrow(catComb)){
#	i0 = eCo[which(eCo$category1==catComb[i,1] & eCo$category2==catComb[i,2]),]
#	kW = kwAllPairsNemenyiTest(i0$ratio_in_rep~as.factor(i0$c1_is), dist="Chisquare", p.adjust.method="none")
#	kwSta = as.data.frame(kW$statistic)
#	kwPva = as.data.frame(fdrBH(kW$p.value))
#	kwS = kwP = c();for(j in 1:ncol(kwSta)){
#		kwS = c(kwS,kwSta[,j])
#		kwP = c(kwP,kwPva[,j])
#	};kwTab = data.frame(
#		"interaction1"=rep(row.names(kwSta),ncol(kwSta)),
#		"interaction2"=rep(colnames(kwSta),each=nrow(kwSta)),
#		"chi.sq"=kwS, "adj.p"=kwP)
#	kwTab = kwTab[!is.na(kwTab$adj.p),]
#	write.csv(kwTab,paste0(pOut,nAm,"-kwPairs_",catComb[i,1],"_",catComb[i,2],".csv"), quote=F, row.names=F)

## grouped boxplot
#	for(i2 in 1:length(tY)){
#		i0$c1_is[which(i0$c1_is==tY[i2])] = tY0[i2]
#	}
#	pdf(paste0(pOut,nAm,"-kwPairs_",catComb[i,1],"_",catComb[i,2],".pdf"), width=7,height=9)
#	par(mar=c(7,4.5,0,2)+.1, xpd=T)
#	boxplot(i0$ratio_in_rep~i0$c1_is,
	#boxplot(i0$ratio_in_rep~gsub("/"," / ",gsub("_"," ",gsub("_c2","",gsub("_of_c2","",i0$c1_is)))),
#	ylim=c(0,1+nrow(kwTab)/5), col="#FFFFFFFF", xlab="", ylab=paste0("ratio of likeliness in ",argv[4]," replicates"), las=2, pch=4, yaxt="n", cex.lab=2, cex.axis=2)
#	axis(2,at=seq(0,1,.2),labels=seq(0,1,.2))
#	if(nrow(kwTab)>0){
#		lB = tY[order(tY)]
#		lBplt = kwTab
#		for(i2 in 1:length(lB)){for(i1 in 1:2){lBplt[which(lBplt[,i1]==lB[i2]),i1]=i2}}
#		lBplt$chi.sq = round(lBplt$chi.sq,2)
#		lBplt$adj.p = ifelse(lBplt$adj.p<.001,"<<0.01",round(lBplt$adj.p,3))
#		cOls = ifelse(lBplt$adj.p>.1,"#00000022","#000000ff")
#		segments(x0=as.numeric(lBplt[,1]),x1=as.numeric(lBplt[,2]),y0=(1:nrow(lBplt))/5+1, col=cOls)
#		text(x = (as.numeric(lBplt[,1])+as.numeric(lBplt[,2]))/2-.5, y = (1:nrow(lBplt))/5+1.1, labels=paste("X =",lBplt[,3],"; adj-p =",lBplt[,4]), xpd=T, col=cOls)
#	}
#	invisible(dev.off())
#}

