#!/bin/bash
# author: ph-u
# script: analysisPlots.sh
# desc: plot + ecology analysis
# in: bash analysisPlots.sh [path/2/data]
# out: NA (see respective scripts)
# arg: 1
# date: 20220111

pT=$1
for i in `ls ${pT}/*-log.csv | rev | cut -f 2 -d "-" | cut -f 1 -d "/" | rev`;do
	echo -e ${i}
	Rscript interactionTypes.r ${pT}/ ${i} 2> /dev/null &
	Rscript SimDataPlot.r ${pT}/ ${i}
done
exit
