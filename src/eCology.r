#!/bin/env Rscript
# author: ph-u
# script: eCology.r
# desc: extract ecology from simulation replicates
# in: Rscript eCology.r [path/to/data/] [time_series_basename] [LVC/gLV] [replicates]
# out: data/*-{eco,tsMatch,rkValue,filter}.csv, result/*-ts{Data,AllRep}.pdf
# arg: 4
# date: 20220508 (supersede interactionTypes.r), 20220820 (CSD3 adaptation)

#SBATCH -J ecolAna
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p skylake-himem

##### env #####
argv=(commandArgs(T))
source("src.r")
library(deSolve)
percCFU = 42; # population threshold to distinguish between CFU and % type data
#library(PMCMRplus) # v1.9.3
#source("fdrBH.r")
pT = argv[1]; nAm = argv[2]; tYpe=argv[3]; pOut = gsub("data","result",pT) # path links
t0 = read.csv(paste0(pT,nAm,"-log.csv"), header=T) # processed time-series data
pR = read.csv(paste0(pT,nAm,"-pri.csv"), header=T) # prior list
sD = read.csv(paste0(pT,nAm,"-seed.csv"), header=T) # seed list
rP = rK = vector(mode="list", length=as.numeric(argv[4])) # bin for simulated coefficients
for(i in 1:as.numeric(argv[4])){
	rP[[i]] = read.csv(paste0(pT,nAm,"-",i,"-sam.csv"), header=T)
}
dMaxSim = nrow(rP[[1]])*as.numeric(argv[4]) # ref max simulation
n=colnames(t0)[-1] # list of categories name
if(any(t0[,-1]>percCFU)){yLab="Presence in pwCF [%]"}else{yLab="log_e(y+1) [CFU/mL]"} # y-axis label
oDe = ifelse(tYpe=="LVC","c","g") # equation type
x0 = rep(0,ncol(t0)-1) # initiate population vector
for(i in 2:ncol(t0)){x0[i-1] = median(t0[which(t0[,1]==min(t0[,1])),i])} # set populations (20220411)
if(R.Version()$major>=4){ # set plot colours
	cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)
	cBl = palette.colors(palette = "Okabe-Ito", alpha=.1, recycle = T)
}else{
	cBp = c("#000000FF","#E69F00FF","#56B4E9FF","#009E73FF","#F0E442FF","#0072B2FF","#D55E00FF","#CC79A7FF","#999999FF")
	cBp = rep(cBp,ceiling(length(n)/length(cBp)))[1:length(n)]
	cBl = paste0(substr(cBp,1,7),"33", sep="")
}

##### categorize ecology ##### (from interactionTypes.r)
#tY0 = c("mutu","cHos","prey","comm","neut","hmED","pred","hmIG","comp")
tY = c("mutualism","commensal_Host_of_c2","prey/host_of_c2","commensal_of_c2","neutral/no_interaction","harmed_by_c2","predator/parasite_of_c2","harming_c2","competition")
sEr = function(i1,i2,cT=tY,tP=tYpe){
        rEf = data.frame("P2p"=rep(1:-1,each=3),"p2P"=rep(1:-1,3),"P_is_the"=cT, stringsAsFactors=F)
        if(tP=="LVC"){i1 = -i1;i2 = -i2}
        i1 = ifelse(i1==0,0,ifelse(i1>0,1,-1))
        i2 = ifelse(i2==0,0,ifelse(i2>0,1,-1))
        iT = rep(NA,length(i1))
        for(i in 1:length(iT)){iT[i] = rEf$P_is_the[which(rEf$P2p==i1[i] & rEf$p2P==i2[i])]}
        return(iT)
}

##### get simulation match range ##### (from SimDataPlot.r)
acRatio = .25
tUq = unique(t0[,1])
tUq = tUq[order(tUq)[-1]] ## accending order safety net
dMin = dMax = dRec = as.data.frame(matrix(0,nr=length(tUq),nc=ncol(t0)))
colnames(dMin) = colnames(dMax) = colnames(dRec) = colnames(t0)
dMin[,1] = dMax[,1] = dRec[,1] = tUq
for(i in 1:length(tUq)){
        d = t0[which(t0[,1]==tUq[i]),]
        for(j in 2:ncol(t0)){
                d0 = range(d[,j])
                dMin[i,j] = max(0, d0[1]-ifelse(length(t0[,1])==length(unique(t0[,1])),acRatio*100,diff(d0)*1.5)) # boxplot outlier definition
                dMax[i,j] = min(100, d0[2]+ifelse(length(t0[,1])==length(unique(t0[,1])),acRatio*100,diff(d0)*1.5))
}}

##### plot legend format ##### https://stackoverflow.com/questions/39552682/base-r-horizontal-legend-with-multiple-rows
nDim = 3;nDim = c(nDim, ceiling((ncol(t0)-1)/nDim))
legMx = matrix(1:prod(nDim), nrow=nDim[1], ncol=nDim[2], byrow=F)
legBd = rep("#000000ff", prod(nDim))
legBd[legMx>(ncol(t0)-1)] = legMx[legMx>(ncol(t0)-1)] = NA

##### set parameter value bin #####
paraBin = as.data.frame(matrix(nr=0,nc=ncol(rP[[1]])+1))
colnames(paraBin) = c("replicate",colnames(rP[[1]]))

##### plot Time-series #####
fCap = function(x){paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")}
colnames(t0) = fCap(colnames(t0))
for(i in c(1:2)){
	z = ifelse(i>1,"AllRep","Data")
	pdf(paste0(pOut,nAm,"-ts",z,".pdf"), width=16)
	par(mar=c(14,5,1,3)+.1, xpd=T) #par(mar=c(5,5,1,12)+.1, xpd=T)
	matplot(t0[,1],t0[,-1], type="p", pch=(1:(ncol(t0)-1))%%25, cex=2, col=cBp,
		xlab=paste0(gsub("_"," (",colnames(t0)[1]),ifelse(length(grep("_",colnames(t0)[1]))>0,")","")),
		ylab=yLab, cex.axis=2, cex.lab=2)
	if(i==1){
		legend("bottom", inset=c(0,-.75), legend = colnames(t0)[-1][legMx], title="Taxonomic Category", border=NA, xpd=T, cex=2, ncol=nDim[2], pch = c((1:(ncol(t0)-1))%%25,NA), lty=c((1:(ncol(t0)-1))%%5+1,NA), lwd=2, col = cBp)
		invisible(dev.off())}
#legend("topright", inset=c(-.19,0), legend = colnames(t0)[-1], pch = (1:(ncol(t0)-1))%%25, lty=(1:(ncol(t0)-1))%%5+1, lwd=2, col = cBp)
}

for(i1 in 1:length(rP)){
	set.seed(sD$seed[i1])
	for(i in 1:nrow(rP[[i1]])){ tK = 0
		a0 = solveLV(x0, as.numeric(rP[[i1]][i,]), range(t0[,1]), oDe)
		a1 = (dMin[,-1]<=a0[which(a0[,1] %in% dMin[,1]),-1]) & (a0[which(a0[,1] %in% dMax[,1]),-1]<=dMax[,-1]); a1[is.na(a1)] = 0 # Simulation-data match count (20220822)
## if CFU, simulation-data match threshold = half
		if(all(colSums(a1)>=(nrow(a1)*ifelse(any(t0[,-1]>percCFU),1,.5)))){tK = 1
			if(nrow(a1)>2 & any(t0[,-1]>percCFU)){ for(i0 in 1:(nrow(a1)-1)){
				if(any(colSums(a1[i0:(i0+1),-1])<2)){tK = 0;break}
			}}}
		if(tK>0){
			paraBin[nrow(paraBin)+1,] = c(i1,as.numeric(rP[[i1]][i,]))
			matplot(a0[,1],a0[,-1], type="l", add=T, lty=(1:(ncol(t0)-1))%%5+1, col=cBl)
		}
		dRec[,-1] = dRec[,-1] + a1
        }#;rP[[i1]] = paraBin[which(paraBin$replicate==i1),-1]
#	rK[[i1]] = rP[[i1]][,which(pR[,2] == "r" | pR[,2] == "k")]
#	rP[[i1]] = rP[[i1]][,which(pR[,2] != "r" & pR[,2] != "k")]
};rm(i,i1)
legend("bottom", inset=c(0,-.75), legend = colnames(t0)[-1][legMx], title=paste("Taxonomic Category -",nrow(paraBin),"simulation set(s)"), border=NA, xpd=T, cex=2, ncol=nDim[2], pch = c((1:(ncol(t0)-1))%%25,NA), lty=c((1:(ncol(t0)-1))%%5+1,NA), lwd=2, col = cBp)
invisible(dev.off())

write.csv(paraBin,paste0(pT,nAm,"-filter.csv"), quote=F, row.names=F)
#pRrk = pR[which(pR[,2] == "r" | pR[,2] == "k"),]
#pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]
dRec[,-1] = dRec[,-1]/dMaxSim # data-matching ratio (from SimDataPlot.r)
write.csv(dRec,paste0(pT,nAm,"-tsMatch.csv"), quote=F, row.names=F)

##### ex: growth rate, carrying capacity #####
rK = which(pR[,2] %in% c("r","k"))
rkVal = paraBin[,c(1,1+rK)]
colnames(rkVal)[-1] = paste(pR[rK,1],pR[rK,2],sep=".")
#if(nrow(rK[[1]])>0){
#	rkVal = rK[[1]]
#	rkVal$replicate = 1
#}else{
#	rkVal = as.data.frame(matrix(0,nr=0,nc=nrow(pR)+1))
#	colnames(rkVal)[ncol(rkVal)]="replicate"
#}
#for(i in 2:length(rK)){if(nrow(rK[[i]])>0){
#	tMp = rK[[i]];tMp$replicate = i
#	rkVal = rbind(rkVal,tMp)
#}};colnames(rkVal)[-ncol(rkVal)] = paste(pRrk[,1],pRrk[,2],sep=".")
write.csv(rkVal,paste0(pT,nAm,"-rkValue.csv"), quote=F, row.names=F)

##### interaction matrix ##### (from interactionTypes.r)
n = as.character(unique(pR[,1]))
iNt = paraBin[,-c(1,1+rK)]
if(tYpe=="LVC"){
        a = as.data.frame(matrix(1, nr=length(n), nc=length(n)), stringsAsFactors=F)
        i0=1
        for(i in 1:nrow(a)){ for(j in 1:ncol(a)){if(i!=j){
                a[i,j] = colnames(iNt)[i0]
                i0=i0+1
}}}}else{
        a = as.data.frame(t(matrix(colnames(iNt), nr=length(n), nc=length(n))), stringsAsFactors=F)
};row.names(a) = colnames(a) = n

##### category pairwise combinations #####
catComb = data.frame(c1=rep(n,each=length(n)),c2=n,c3=NA, stringsAsFactors=F)
for(i in 1:nrow(catComb)){catComb$c3[i] = paste(catComb[i,-3][order(catComb[i,-3])],collapse=".")}
catComb = catComb[!duplicated(catComb$c3),-3]
if(tYpe=="LVC"){catComb = catComb[which(catComb$c1!=catComb$c2),]}
#catComb = as.data.frame(matrix(NA,nr=choose(length(n),2)+ifelse(tYpe=="LVC",0,length(n)), nc=2))
#nR=i=j=1;repeat{
#	if(ifelse(tYpe=="LVC",i<j,i<=j)){catComb[nR,] = c(n[i],n[j]);nR = nR+1}
#	j = j+1
#	if(j>length(n)){j=1;i = i+1}
#	if(nR>nrow(catComb)){break}
#};rm(nR,i,j)

##### interaction categorization result collector #####
rEp = unique(paraBin$replicate)
eCo = as.data.frame(matrix(0,nr=nrow(catComb)*length(rEp)*length(tY), nc=7), stringsAsFactors=F)
#eCo = as.data.frame(matrix(0,nr=(choose(length(n),2)+ifelse(tYpe=="LVC",0,length(n)))*length(rP)*length(tY),nc=7))
colnames(eCo) = c(paste0("category",1:2),"replicate","c1_is","count","fit_sim","ratio_in_rep")
eCo$replicate = rep(rEp,each=nrow(eCo)/length(rEp))#rep(1:length(rP),each=nrow(eCo)/length(rP))
eCo$c1_is = rep(tY,nrow(eCo)/length(tY))
for(i in 1:2){eCo[,i] = rep(rep(catComb[,i],each=length(tY)),length(rEp))}

##### map interactions #####
for(i in 1:length(rEp)){ #for(i in 1:length(rP)){
	a0 = paraBin[which(paraBin$replicate==rEp[i]),-c(1,1+rK)]
	eCo$fit_sim[which(eCo$replicate==rEp[i])] = nrow(a0) #nrow(rP[[i]])
	for(c2 in 1:length(n)){ for(c1 in 1:length(n)){
		if(ifelse(tYpe=="LVC",c2>c1,c2>=c1)){
			cAt = table(sEr(a0[,a[c2,c1]],a0[,a[c1,c2]])) # a[row,col] - "col" (P) affect population of "row" (p)
			for(i0 in 1:length(cAt)){
				eCo$count[which(eCo$category1==colnames(a)[c1] & eCo$category2==colnames(a)[c2] & eCo$replicate==rEp[i] & eCo$c1_is==names(cAt)[i0])] = cAt[i0]
}
		}
	}}
};rm(i)
eCo$ratio_in_rep = eCo$count/eCo$fit_sim # relationship ratio in the top 100 best-fit after double biological simulation-data reality check
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

