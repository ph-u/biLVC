#!/bin/env Rscript
# author: ph-u
# script: SimDataPlot.r
# desc: estimated parameter changes through time
# in: p="path/2/data";for i in `ls ${p}/*-log.csv | cut -f 1 -d "-" | rev | cut -f 1 -d "/" | rev`;do Rscript SimDataPlot.r ${p} ${i};done
# out: result/*-{ts}.pdf
# arg: 2
# date: 20220109

argv=(commandArgs(T))
source("src.r")
library(deSolve)
pT = argv[1]
nAm = argv[2]
#nAm = "2_PAO1_SA25923"
p = read.csv(paste0(pT,nAm,"-sam.csv"), header=T, stringsAsFactors=F)
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T)

x0 = rep(0,ncol(t0)-1)
for(i in 2:ncol(t0)){
	x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])
}

##### colour #####
cBp = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#e79f00", "#9ad0f3", "#F0E442", "#999999", "#cccccc", "#6633ff", "#00FFCC", "#0066cc", "#000000")
cBl = c("#E69F0033", "#56B4E933", "#009E7333", "#0072B233", "#D55E0033", "#CC79A733", "#e79f0033", "#9ad0f333", "#F0E44233", "#99999933", "#cccccc33", "#6633ff33", "#00FFCC33", "#0066cc33", "#00000033")
#n = c("PAO1","SA25923","SC5314")
n=colnames(t0)[-1]
cOl = data.frame(id=n,cPt=cBp[1:length(n)],cLn=cBl[1:length(n)])
ptCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),2]
lnCol = cOl[which(cOl[,1] %in% colnames(t0)[-1]),3]

##### plot time-series #####
pdf(paste0(pT,"../result/",nAm,"-ts.pdf"), width=14)
par(mar=c(5,4,0,0)+.1, xpd=F)
matplot(t0[,1],t0[,-1], type="p", pch=1:(ncol(t0)-1), cex=.3, col=ptCol, xlab=colnames(t0)[1], ylab="log(y+1) [CFU/mL]")
legend("bottomright", inset=c(0,0), legend = colnames(t0)[-1], pch = rep(16,ncol(t0)-1), col = ptCol)

for(i in 1:nrow(p)){
	a0 = solveLVC(x0, as.numeric(p[i,]), range(t0[,1]))
	matplot(a0[,1],a0[,-1], type="l", add=T, col=lnCol)
}
invisible(dev.off())
