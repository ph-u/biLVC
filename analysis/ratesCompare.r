#!/bin/env Rscript
# author: ph-u
# script: ratesCompare.r
# desc: growth rate & carrying capacity comparisons
# in: none
# out: none
# arg: NA
# date: 20220110

psc = read.csv("../data/3_PAO1_SA25923_SC5314-sam.csv", header=T)
ps = read.csv("../data/2_PAO1_SA25923-sam.csv", header=T)
pc = read.csv("../data/2_PAO1_SC5314-sam.csv", header=T)
sc = read.csv("../data/2_SA25923_SC5314-sam.csv", header=T)
rEf2 = rep(c("r","k","a"),2)
rEf3 = rep(c("r","k","a","a"),3)
n = c("PAO1","SA25923","SC5314")

##### growth rates #####
psc_r = psc[,which(rEf3=="r")]
ps_r = ps[,which(rEf2=="r")]
pc_r = pc[,which(rEf2=="r")]
sc_r = sc[,which(rEf2=="r")]
colnames(psc_r) = n
colnames(ps_r) = n[-3]
colnames(pc_r) = n[-2]
colnames(sc_r) = n[-1]

i=1
pdf(paste0("../result/growthRate_",n[i],".pdf"))
boxplot(c(psc_r[,n[i]],ps_r[,n[i]],pc_r[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="growth rate [CFU/mL/hr]", xlab="Co-culture experiment", ylim=c(-.3,2.1))
p = t.test(psc_r[,n[i]],ps_r[,n[i]]) 
segments(2,-.1,3,-.1);text(2.5,-.15,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_r[,n[i]],pc_r[,n[i]]) 
segments(1,-.3,3,-.3);text(2,-.35,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(pc_r[,n[i]],ps_r[,n[i]]) 
segments(1,-.2,2,-.2);text(1.5,-.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

i=2
pdf(paste0("../result/growthRate_",n[i],".pdf"))
boxplot(c(psc_r[,n[i]],ps_r[,n[i]],sc_r[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="growth rate [CFU/mL/hr]", xlab="Co-culture experiment", ylim=c(-.3,1.4))
p = t.test(psc_r[,n[i]],ps_r[,n[i]]) 
segments(2,-.1,3,-.1);text(2.5,-.13,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_r[,n[i]],sc_r[,n[i]]) 
segments(1,-.3,3,-.3);text(2,-.33,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(sc_r[,n[i]],ps_r[,n[i]]) 
segments(1,-.2,2,-.2);text(1.5,-.23,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

i=3
pdf(paste0("../result/growthRate_",n[i],".pdf"))
boxplot(c(psc_r[,n[i]],pc_r[,n[i]],sc_r[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="growth rate [CFU/mL/hr]", xlab="Co-culture experiment", ylim=c(-.2,.4))
p = t.test(psc_r[,n[i]],pc_r[,n[i]]) 
segments(2,-.05,3,-.05);text(2.5,-.07,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_r[,n[i]],sc_r[,n[i]]) 
segments(1,-.1,3,-.1);text(2,-.12,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(pc_r[,n[i]],sc_r[,n[i]]) 
segments(1,-.15,2,-.15);text(2.5,-.17,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

##### carrying capacities #####
psc_k = psc[,which(rEf3=="k")]
ps_k = ps[,which(rEf2=="k")]
pc_k = pc[,which(rEf2=="k")]
sc_k = sc[,which(rEf2=="k")]
colnames(psc_k) = n
colnames(ps_k) = n[-3]
colnames(pc_k) = n[-2]
colnames(sc_k) = n[-1]

i=1;l = 13
pdf(paste0("../result/carryingCapacity_",n[i],".pdf"))
boxplot(c(psc_k[,n[i]],ps_k[,n[i]],pc_k[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="carrying capacity [CFU/mL]", xlab="Co-culture experiment", ylim=c(l-1,22))
p = t.test(psc_k[,n[i]],ps_k[,n[i]]) 
segments(2,l+1,3,l+1);text(2.5,l+.75,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_k[,n[i]],pc_k[,n[i]]) 
segments(1,l,3,l);text(2,l-.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(pc_k[,n[i]],ps_k[,n[i]]) 
segments(1,l+.5,2,l+.5);text(1.5,l+.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

i=2;l=10
pdf(paste0("../result/carryingCapacity_",n[i],".pdf"))
boxplot(c(psc_k[,n[i]],ps_k[,n[i]],sc_k[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="carrying capacity [CFU/mL/hr]", xlab="Co-culture experiment", ylim=c(8,23))
p = t.test(psc_k[,n[i]],ps_k[,n[i]]) 
segments(2,l+1,3,l+1);text(2.5,l+.75,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_k[,n[i]],sc_k[,n[i]]) 
segments(1,l,3,l);text(2,l-.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(sc_k[,n[i]],ps_k[,n[i]]) 
segments(1,l+.5,2,l+.5);text(1.5,l+.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

i=3;l=0
pdf(paste0("../result/carryingCapacity_",n[i],".pdf"))
boxplot(c(psc_k[,n[i]],pc_k[,n[i]],sc_k[,n[i]])~rep(c("PSC","PS","PC"), each=nrow(psc)), ylab="carrying capacity [CFU/mL/hr]", xlab="Co-culture experiment", ylim=c(l,17))
p = t.test(psc_k[,n[i]],pc_k[,n[i]]) 
segments(2,l+1,3,l+1);text(2.5,l+.75,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(psc_k[,n[i]],sc_k[,n[i]]) 
segments(1,l,3,l);text(2,l-.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
p = t.test(pc_k[,n[i]],sc_k[,n[i]]) 
segments(1,l+.5,2,l+.5);text(1.5,l+.25,paste("T =",round(p$statistic,2),"; adj.p =",round(p$p.value*length(n),4)))
invisible(dev.off())

