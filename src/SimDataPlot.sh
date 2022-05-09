#!/bin/bash
# author: ph-u
# script: SimDataPlot.sh
# desc: plot selected simulation replicate result to data
# in: bash analysisPlots.sh [/full/path/2/*-rec.txt]
# out: NA (see respective scripts)
# arg: 1 (req: *-{sam,pri,log}.csv, *-rec.txt)
# date: 20220509 (segregate from analysisPlots.sh)

[[ -z $1 ]] && grep -e "^\# desc\|^\# in" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit 
pT=`dirname $1`
bNam=`echo -e $1 | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "-"`
tYpe=`echo -e $1 | rev | cut -f 1 -d "_" | rev | cut -f 1 -d "-"`
rEp=`echo -e $1 | rev | cut -f 2 -d "-" | rev`
nAcc=`grep -e "accepted" $1 | tail -n 1 | cut -f 2 -d ":" | cut -f 2 -d " "`
Rscript SimDataPlot.r ${pT}/ ${bNam} ${rEp} ${tYpe} ${nAcc}
exit
