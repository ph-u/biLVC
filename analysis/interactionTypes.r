#!/bin/env Rscript
# author: ph-u
# script: interactionTypes.r
# desc: interaction types through time
# in: p="path/2/data";for i in `ls ${p}/*-log.csv | cut -f 1 -d "-" | rev | cut -f 1 -d "/" | rev`;do Rscript interactionTypes.r ${p} ${i};done
# out: result/*-eco*.csv
# arg: 2
# date: 20220109

argv=(commandArgs(T))
pT = argv[1]
nAm = argv[2]
p0 = read.csv(paste0(pT,nAm,"-sam.csv"), header=T, stringsAsFactors=F)
pR = read.csv(paste0(pT,nAm,"-pri.csv"), header=T, stringsAsFactors=F)

p0 = p0[,which(pR[,2] != "r" & pR[,2] != "k")]
pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]

##### single ecological relationship ##### 2022{0109,0111}
tY = c("competition","mutualism","neutral","predator/parasite","prey/host","commensal","commensal_Host","harmed_by","harming_other")
sEr = function(i1,i2,cT=tY){ # + harm ; - help ; compare i1 (1->2) to i2 (2->1)
	iT = ifelse(i1<0, ifelse(i2<0,tY[2], ifelse(i2>0,tY[5],tY[7])),
	ifelse(i1>0, ifelse(i2>0,tY[1],ifelse(i2<0,tY[4],tY[9])),
	ifelse(i2>0,tY[8],ifelse(i2<0,tY[6],tY[3]))))
	return(iT)
}

##### map columns into interaction matrix #####
n = unique(pR[,1])
a = as.data.frame(matrix(1, nr=length(n), nc=length(n)))
row.names(a) = colnames(a) = n
i0=1
for(i in 1:nrow(a)){ for(j in 1:ncol(a)){if(i!=j){
	a[i,j] = colnames(p0)[i0]
	i0=i0+1
}}}

##### proportion of deduced ecological interaction types #####
w = as.data.frame(matrix(0, nr=ncol(p0)/2, nc=length(tY)+4))
colnames(w) = c("subject","other",tY,"Chi-sq","p-val")
i0=1; for(i in 2:length(n)){for(j in 1:(length(n)-1)){if(i>j){
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
				e[j0,] = c(tY[c(i1,j1)],b0$statistic,b0$p.value*nrow(e))
				j0=j0+1}}}
			write.csv(e, paste0(pT,"../result/",nAm,"-eco_",n[i],"_",n[j],".csv"), quote=F, row.names=F)
	}}else{w[i0,(ncol(w)-1):ncol(w)] = rep(NA,2)};i0=i0+1}}}
write.csv(w, paste0(pT,"../result/",nAm,"-eco.csv"), quote=F, row.names=F)
