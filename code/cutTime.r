#!/bin/env Rscript
# author: ph-u
# script: cutTime.r
# desc: cut time-series into sections
# in: Rscript cutTime.r [basename] [cut-time...]
# out: raw/[basename]_[start][end].csv
# arg: >=2 (2nd onwards are cut-time(s))
# date: 20220112

argv=(commandArgs(T))
nAm = argv[1]
z = read.csv(paste0("../raw/",nAm,".csv"), header=T, stringsAsFactors=F)
cUt = c(range(z[,1]),as.numeric(argv[-1]))
cUt = cUt[order(cUt)]
for(i in 2:length(cUt)){
	o = z[which(z[,1]>=cUt[i-1] & z[,1]<=cUt[i]),]
#	if(min(o[,1])>0){o[,1] = o[,1]-min(o[,1])}
	s = paste0(paste0(rep(0,nchar(cUt[length(cUt)])-nchar(cUt[i-1])),collapse=""),cUt[i-1])
	e = paste0(paste0(rep(0,nchar(cUt[length(cUt)])-nchar(cUt[i])),collapse=""),cUt[i])
	write.csv(o,paste0("../raw/",nAm,"_",s,e,".csv"), quote=F, row.names=F)
}
