#!/bin/env Rscript
# author: ph-u
# script: src.r
# desc: functions for bayesInfer.r
# in: source("src.r")
# out: NA
# arg: 0
# date: 20220108

##### parameter arrangement ##### 20220107
pAra = function(p, n){
	t0 = matrix(p,nc=n)
	t1 = as.matrix(t(t0[3:nrow(t0),]))
	a0 = matrix(1,nr=n,nc=n)
	if(n>2){
		for(i in 1:nrow(t1)){for(j in 1:ncol(t1)){if(j<i){a0[i,j] = t1[i,j]}else{a0[i,j+1] = t1[i,j]}}}
	}else{a0[1,2] = t1[1,1];a0[2,1] = t1[1,2]}
	return(list(r=t0[1,], k=t0[2,], a=a0))
}

##### Lotka-Volterra Competition (LVC) Model ##### 20220107
LVc = function(t, x, pAr){with (as.list(c(x, pAra(pAr,length(x)))), {
	iNt = 0
	for(i in 1:length(x)){iNt = iNt + x[i] * a[,i]}
	return(list(r * x * (1 - iNt / k)))
	})}

##### Simulation ##### 20220107
solveLVC = function(u, pAr, tspan, mEthod="rk4"){return(as.data.frame(ode(y=u, times=seq(tspan[1], tspan[2], by=1), func=LVc, parms=pAr, method=mEthod)))}

##### Objective function to minimise, >=2 components ##### 20220107
oBj = function(pAr, df, mEthod="rk4"){
	tspan = range(df[,1])
	u0 = df[which(df[,1]==tspan[1]),-1]
	u = rep(0,ncol(u0))
	for(i in 1:length(u)){u[i]=median(u0[,i])}
	s = solveLVC(u, pAr, tspan, mEthod)
	colnames(s) = colnames(df)
	cOst = modCost(obs = df[,c(1,2)], model = s[,c(1,2)], x=colnames(df)[1])
	for(i in 3:ncol(df)){cOst = modCost(obs = df[,c(1,i)], model = s[,c(1,i)], cost=cOst, x=colnames(df)[1])}
	return(cOst)
}

##### Compare Integration Methods ##### 20220107
bestMet = function(pAr, df){
	#fIts = data.frame(method=c("euler", "rk4", "lsoda"), fitness=rep(0,3), stringsAsFactors=F)
	fIts = data.frame(method=c("rk4", "euler"), fitness=rep(0,2), stringsAsFactors=F)
	for(i in 1:nrow(fIts)){fIts[i,2] = oBj(pAr, df, fIts[i,1])$model}
	return(fIts[which(fIts[,2]==min(fIts[,2])),1][1])
}

##### Model residuals for MCMC ##### 20220108
mCres = function(pAr, df){return(oBj(pAr, df, bestMet(pAr, df))$model)}

##### rolling growth rate ##### 20220101
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
	return(c(rEs[which(rEs[,3]==max(rEs[,3])),1], p1+mean(p1)*c(-1,1)))
}
