#!/bin/bash
# author: ph-u
# script: analysisPlots.sh
# desc: plot + ecology analysis
# in: bash analysisPlots.sh [path/2/data]
# out: NA (see respective scripts)
# arg: 1
# date: 20220111

[[ -z $1 ]] && grep -e "^\# desc\|^\# in" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit
pT=$1
for i in `ls ${pT}/*-log.csv | rev | cut -f 2 -d "-" | cut -f 1 -d "/" | rev`;do
#	echo -e ${i}
	tYpe=`echo ${i} | rev | cut -f 1 -d "_" | rev`
	nAcc=`grep -e "accepted" ${pT}/${i}-rec.txt | cut -f 2 -d ":" | cut -f 2 -d " "`
	Rscript interactionTypes.r ${pT}/ ${i} ${tYpe} ${nAcc} 2> /dev/null &
	Rscript SimDataPlot.r ${pT}/ ${i} ${tYpe} ${nAcc}
done
exit
