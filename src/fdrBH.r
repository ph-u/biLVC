#!/bin/env Rscript
# author: ph-u
# script: fdrBH.r
# desc: Benjamini-Hochberg False Discovery Rate correction wrapper
# in: source("fdrBH.r")
# out: NA
# arg: 0
# date: 20220610

fdrBH = function(pMtx){
	z = matrix(p.adjust(pMtx,method="BH"),nc=ncol(pMtx))
	colnames(z) = colnames(pMtx)
	row.names(z) = row.names(pMtx)
	return(z)
}
