#!/bin/env Rscript
# author: ph-u
# script: ecoBars.r
# desc: barplot of ecological distribution
# in: Rscript ecoBars.r [path/2/data] [sampleFreq_basename] [replicate] [simulation_output_size] [log10?]
# out: ?
# arg: 0
# date: 20230101

##### in #####
argv = (commandArgs(T))
ptIN = paste0(argv[1],"/"); ptOT = sub("data","result",ptIN)
nRep = as.numeric(argv[3]); simO = as.numeric(argv[4])
mTx = .3
pFx = strsplit(argv[2],"_")[[1]][1] # freq file prefix = ecoPlots'
sEq = c(2,9,3,4,5,6,7,1,8)
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[sEq] # >= R v4.0.0
sFreq = read.csv(paste0(ptIN,list.files(ptIN,argv[2])), header=T) # sample frequency file must be unique

##### f: capitalise first letter #####
capFirst = function(x){return(paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))))}

##### f: plot legend #####
legPlot = function(x,nDim=3){
        nDim = c(nDim, ceiling(length(x)/nDim))
        legMx = matrix(1:prod(nDim), nrow=nDim[1], ncol=nDim[2], byrow=F)
        legBd = rep("#000000ff", prod(nDim))
        legBd[legMx>length(x)] = legMx[legMx>length(x)] = NA
        return(list(legMx,nDim[2]))
}

##### prescription popularity #####
cat("grouping sample sizes plot: ",date(),"\n")
pdf(paste0(ptOT,pFx,"_samFreq.pdf"), width=14)
par(mar=c(5,4.5,1.2,1)+.1, mfrow=c(2,1),xpd=F)
if(any(c("Y","y","YES","yes","Yes","log","LOG","Log","log10","LOG10","Log10") %in% argv[5])){z = log10(sFreq[,-1]); z0 = "log10(Sample)"}else{z = sFreq[,-1]; z0 = "Sample"}
matplot(sFreq[,1],z,type="b",pch=1:(ncol(sFreq)-1),lty=1,lwd=3,col=cBp, xlab=capFirst(colnames(sFreq)[1]), ylab=z0, cex.axis = 2.1, cex.lab=2.1)
rm(z,z0)
plot.new();lPt = legPlot(colnames(sFreq)[-1], 2)
legend("top",legend=capFirst(colnames(sFreq)[-1])[lPt[[1]]], border=NA, ncol=lPt[[2]], lty=c(rep(1,length(1:(ncol(sFreq)-1))),NA), title="Groups", col=cBp[1:(ncol(sFreq)-1)], lwd=5, cex=2.1, pch=c(1:(ncol(sFreq)-1),NA))
invisible(dev.off())

##### files #####
cat("Import file list: ",date(),"\n")
f = list.files(ptIN,"-eco")
f0 = unlist(strsplit(f,"_"))
fNam = as.data.frame(matrix(nr=length(f),nc=2))
# file format: [prefix]_[group]_[startEnd]_[ecoEquation]-eco.csv
for(i in 1:nrow(fNam)){fNam[i,] = f0[(2:3)+length(f0)/length(f)*(i-1)]};rm(i,f0)
fNam$start = substr(fNam[,2],1,2); fNam$end = substr(fNam[,2],3,4)
f0 = unique(fNam[,1])
## time-scale standardization
yR = as.numeric(unique(c(fNam$start,fNam$end)))+ifelse(any(c("Year","YEAR") %in% colnames(sFreq)[1]),2000,0)
yR = yR[order(yR)]

##### data bin to summary rolling-time graphs #####
cat("Processing group:\n")
for(i in 1:length(f0)){z = 0; cat(i,"(",date(),"),\n")
        for(i0 in which(fNam[,1]==f0[i])){
                f1 = read.csv(paste0(ptIN,f[i0]), header=T, stringsAsFactors=F)
                t0 = as.numeric(fNam$end[i0])+ifelse(any(c("Year","YEAR") %in% colnames(sFreq)[1]),2000,0)
## set-up record data frame
                if(z==0){
                        if(i==1){
                                sPair = unique(paste0(f1$category1,"_",f1$category2))
                                eCo = unique(f1$c1_is)
                        }
                        ecoTS = vector("list", length(sPair)) # list of pairwise time-series
                        for(i1 in 1:length(ecoTS)){
                                ecoTS[[i1]] = as.data.frame(matrix(0,nr=length(yR)-2, nc=length(eCo)+1))
                                colnames(ecoTS[[i1]]) = c(colnames(sFreq)[1],eCo)
                                ecoTS[[i1]]$year = yR[-c(1:2)]
                        };rm(i1)
                z = 1}
## segregate data (one end time) into dataframes
                for(i1 in 1:length(sPair)){
                        a0 = strsplit(sPair[i1],"_")[[1]]
                        f2 = f1[which(f1$category1==a0[1] & f1$category2==a0[2]),]
                        for(i2 in eCo){ecoTS[[i1]][which(ecoTS[[i1]]$year==t0),which(colnames(ecoTS[[i1]])==i2)] = sum(f2$count[which(f2$c1_is==i2)])} # count one eco relationship
                };rm(i1)
        };rm(i0)
## plot dataframes
        for(i0 in 1:length(sPair)){
                ecoTS[[i0]][,-1] = ecoTS[[i0]][,-1]/(simO*nRep)*100
### stacked barplot prep
                ecoPlot = ecoTS[[i0]][,-1]
                ecoPlot[nrow(ecoPlot)+c(1:2),] = 0
                ecoPlot = as.matrix(t(ecoPlot[c(nrow(ecoPlot)-1:0,1:(nrow(ecoPlot)-2)),]))
                colnames(ecoPlot) = paste0(sFreq[,1],"\n(",sFreq[,which(colnames(sFreq)==f0[i])],")\n(",apply(sFreq[,-1],1,sum),")")

### export
                pdf(paste0(ptOT,pFx,"_",f0[i],"_",sPair[i0],".pdf"), width=14, height=11)
                par(mar=c(5,4.5,2,1)+.1, mfrow=c(2,1), cex.axis=1.4, xpd=T)
                barplot(ecoPlot, ylim=c(0,100), xaxt="n", ylab=paste(sPair[i0], "(%)"), xlab="", col=cBp, border="white", cex.axis=2.1, cex.lab=1.5)
                axis(1, at=-.5+1.2*(1:length(yR)), padj=.7, labels=colnames(ecoPlot))
                mtext(paste0(capFirst(colnames(sFreq)[1])," (Group Sample Size) (Total Sample Size)"),side=1,padj=4.9,cex=2.1)
### legend plot
                plot.new()
                lPt = legPlot(eCo)
                legend("top", inset=c(0,0), legend = capFirst(gsub("_"," ",sub("c2","B",eCo[lPt[[1]]]))), title=paste("Ecological Relationship - 100% =",nRep*simO,"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,length(eCo)), col = cBp)
                invisible(dev.off())
        };rm(i0)
};rm(i);cat("Done",date(),"\n")

