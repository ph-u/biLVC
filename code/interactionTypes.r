#!/bin/env Rscript
# author: ph-u
# script: interactionTypes.r
# desc: interaction types through time
# in: Rscript interactionTypes.r [basename]
# out: stdout text
# arg: 1
# date: 20220109

argv=(commandArgs(T))
nAm = argv[1]
#nAm = "2_PAO1_SA25923"
#nAm = "2_PAO1_SC5314"
#nAm = "2_SA25923_SC5314"
#nAm = "3_PAO1_SA25923_SC5314"
p0 = read.csv(paste0("../data/",nAm,"-sam.csv"), header=T, stringsAsFactors=F)
pR = read.csv(paste0("../data/",nAm,"-pri.csv"), header=T, stringsAsFactors=F)

p0 = p0[,which(pR[,2] != "r" & pR[,2] != "k")]
pR = pR[which(pR[,2] != "r" & pR[,2] != "k"),]

##### single ecological relationship ##### 20220109
sEr = function(i1,i2){ # + harm ; - help ; compare i1 (1->2) to i2 (2->1)
	iT = ifelse(i1<0 & i2<0, "mutualism",
	ifelse(i1>0 & i2>0, "competition",
	ifelse(i1>i2, "predator/parasite", "prey/host")))
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
w = as.data.frame(matrix(0, nr=ncol(p0)/2, nc=8))
colnames(w) = c("subject","other", "competition","mutualism","predator/parasite","prey/host", "Chi-sq","p-val")
i0=1
for(i in 2:length(n)){for(j in 1:(length(n)-1)){if(i>j){
	w[i0,1:2] = c(n[i],n[j])
	b = table(sEr(p0[,a[i,j]],p0[,a[j,i]]))
	w[i0,names(b)] = b
	l = chisq.test(b)
	w[i0,7:8] = c(l$statistic,l$p.value)
	if(length(b)>2 & l$p.value<.05){ ## pairwise Chi-sq tests with bonferroni p-value correction
		e = as.data.frame(matrix(NA, nr=choose(length(b),2), nc=4))
		colnames(e) = c("interaction1","interaction2",rev(rev(colnames(w))[1:2]))
		j0=1
		for(i1 in 1:(length(b)-1)){for(j1 in 2:length(b)){if(i1<j1){
			b0 = chisq.test(b[c(i1,j1)])
			e[j0,] = c(names(b)[c(i1,j1)],b0$statistic,b0$p.value*nrow(e))
			j0=j0+1
		}}}
		write.csv(e, paste0("../result/",nAm,"-eco_",n[i],"_",n[j],".csv"), quote=F, row.names=F)
	}
	i0=i0+1
}}}
write.csv(w, paste0("../result/",nAm,"-eco.csv"), quote=F, row.names=F)
