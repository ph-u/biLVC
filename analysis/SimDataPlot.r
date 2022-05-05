#!/bin/env Rscript
# author: ph-u
# script: SimDataPlot.r
# desc: estimated parameter changes through time
# in: Rscript SimDataPlot.r [path/to/data/] [data_basename] [math-model-type] [num_accepted_simulations]
# out: result/*-{ts}.pdf
# arg: 4
# date: 20220109

argv=(commandArgs(T))
source("../pipeline/src.r")
library(deSolve)
pT = argv[1]; nAm = argv[2]; tYp = argv[3]
#nAm = "2_PAO1_SA25923"
p = read.csv(paste0(pT,nAm,"-sam.csv"), header=T, stringsAsFactors=F)
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T)
nAcc = min(as.numeric(argv[4]), nrow(p))

x0 = rep(0,ncol(t0)-1)
for(i in 2:ncol(t0)){
	x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])
}

##### colour #####
cBp = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#e79f00", "#9ad0f3", "#F0E442", "#999999", "#cccccc", "#6633ff", "#00FFCC", "#0066cc", "#000000")
cBl = c("#E69F0028", "#56B4E944", "#009E7349", "#0072B233", "#D55E0033", "#CC79A733", "#e79f0033", "#9ad0f333", "#F0E44233", "#99999933", "#cccccc33", "#6633ff33", "#00FFCC33", "#0066cc33", "#00000033")
#n = c("PAO1","SA25923","SC5314")
n=colnames(t0)[-1]
cOl = data.frame(id=n,cPt=cBp[1:length(n)],cLn=cBl[1:length(n)])
ptCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),2]
lnCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),3]

##### get simulation match range #####
tUq = unique(t0[,1])
tUq = tUq[order(tUq)[-1]] ## accending order safety net
dMin = dMax = dRec = as.data.frame(matrix(0,nr=length(tUq),nc=ncol(t0)))
colnames(dMin) = colnames(dMax) = colnames(dRec) = colnames(t0)
dMin[,1] = dMax[,1] = dRec[,1] = tUq
for(i in 1:length(tUq)){
	d = t0[which(t0[,1]==tUq[i]),]
	for(j in 2:ncol(t0)){
		d0 = range(d[,j])
		dMin[i,j] = ifelse(length(t0[,1])==length(unique(t0[,1])),min(d0[1]*.9,d0[1]-1),max(d0[1]-diff(d0)/2,0))
		dMax[i,j] = ifelse(length(t0[,1])==length(unique(t0[,1])),max(d0[2]*1.1,d0[2]+1),d0[2]+diff(d0)/2)
}}

##### plot time-series #####
if(any(t0[,-1]>30)){yLab="percentage presence [%]"}else{yLab="log_e(y+1) [CFU/mL]"}
if(argv[3]=="LVC"){oDe="c"}else{oDe="g"}
if(ncol(t0)-1>5){pltAttr = c(T,"topright",-.19,12)}else{pltAttr = c(F,"bottomleft",0,0)}
ptCol0 = rep(ptCol,ceiling((ncol(t0)-1)/length(ptCol)))
lnCol0 = rep(lnCol,ceiling((ncol(t0)-1)/length(lnCol)))

pdf(paste0(pT,"../result/",nAm,"-ts.pdf"), width=14)
par(mar=c(5,5,1,as.numeric(pltAttr[4]))+.1, xpd=as.logical(pltAttr[1]))
matplot(t0[,1],t0[,-1], type="p", pch=(1:(ncol(t0)-1))%%25, cex=1.2, col=ptCol0,
	xlab=paste0(gsub("_"," (",colnames(t0)[1]),ifelse(length(grep("_",colnames(t0)[1]))>0,")","")),
	ylab=yLab, cex.axis=2, cex.lab=2)
legend(pltAttr[2], inset=c(as.numeric(pltAttr[3]),0), legend = colnames(t0)[-1], pch = (1:(ncol(t0)-1))%%25, lty=(1:(ncol(t0)-1))%%5+1, lwd=2, col = ptCol0)

##### plot + simulation percentage match on data #####
i9=nrow(p);for(i in 1:nrow(p)){
	a0 = solveLV(x0, as.numeric(p[i,]), range(t0[,1]), oDe)
	for(i0 in 2:ncol(a0)){a0[,i0] = ifelse(a0[,i0]>150 | a0[,i0]<0,-100,a0[,i0])};rm(i0)
	a0[is.na(a0)] = -100
	if(all(a0!=-100)){
		matplot(a0[,1],a0[,-1], type="l", add=T, col=lnCol0, lty=(1:(ncol(t0)-1))%%5+1)
	}else{i9 = i9-1}
	a0 = a0[which(a0[,1] %in% tUq),]
	for(j in 1:length(tUq)){for(k in 2:ncol(t0)){
		dRec[j,k] = dRec[j,k] + (a0[j,k]>=dMin[j,k] & a0[j,k]<=dMax[j,k])
	}}}
text(fivenum(t0[,1])[2],max(t0[,-1]),paste("Number of simulation(s)\nAccepted:",ifelse(as.numeric(argv[4])>nrow(p), paste0(">",nrow(p)), argv[4]), "; Plotted:",i9), cex=1.2)
invisible(dev.off())

##### simulation percentage match on data #####
dRec[,-1] = dRec[,-1]/nrow(p)
write.csv(dRec,paste0(pT,"../result/",nAm,"-tsMatch.csv"), quote=F, row.names=F)
