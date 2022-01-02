#!/bin/env Rscript
# author: ph-u
# script: paraEst.r
# desc: parameter estimation
# in: Rscript paraEst.r [basename]
# out: data/[basename]-{log,pri}.csv
# arg: 0
# date: 20220101

##### func: rolling growth rate ##### source("func.r")
rolRate = function(df,y,x){
	y = ifelse(is.numeric(y),y,which(colnames(df)==y))
	x = ifelse(is.numeric(x),x,which(colnames(df)==x))
	df = df[,c(x,y)]
	t = unique(df[,x])
	t = t[order(t)] ## accending timeline
	rEs = as.data.frame(matrix(NA, nr=length(t)-1, nc=3))
	for(i in 2:length(t)){
		df0 = df[which(df[,x] %in% c(t[i-1],t[i])),]
## correlation coefficient
		p = summary(lm(df0[,2]~df0[,1]))
		p0 = as.data.frame(p$coefficients)
		rEs[i-1,] = c(p0[2,"Estimate"],p0[2,"Estimate"]*(t[i]-t[i-1]), p$adj.r.squared) # standardise unit time
	}
	p1 = range(rEs[,-3])
	return(c(rEs[which(rEs[,3]==max(rEs[,3])),1], sd(rEs[,-3]), p1+mean(p1)*c(-1,1)))
}

##### import #####
argv=(commandArgs(T))
q = read.csv(paste0("../raw/",argv[1],".csv"), header=T, stringsAsFactors=F)

##### log-transformation #####
q[,-1] = log(q[,-1]+1)
write.csv(q, paste0("../data/",argv[1],"-log.csv"), quote=F, row.names=F)

##### set priors #####
nSp = ncol(q)-1
rEc = as.data.frame(matrix(NA,nr=nSp*(nSp+2),nc=7))
colnames(rEc) = c("id","influencer","par1","par2","type","min","max")
rEc[,1] = rep(colnames(q)[-1], each=nSp+2)
rEc[,2] = rep(c("r","k",colnames(q)[-1]),nSp)
for(i in 1:nrow(rEc)){
	if(rEc[i,2]=="r"){
		p = rolRate(q,rEc[i,1],1)
	}else if(rEc[i,2]=="k"){
		p = range(q[,rEc[i,1]])
		p = c(0,0,p+c(-.2,.2)*mean(p))
	}else{
		p = c(0,1,-2,2)
	}
	rEc[i,3:ncol(rEc)] = c(p[1:2],ifelse(rEc[i,2]=="k","Uniform","Normal"),p[3:4])
}
write.csv(rEc, paste0("../data/",argv[1],"-pri.csv"), quote=F, row.names=F)
