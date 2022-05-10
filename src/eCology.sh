#!/bin/bash
# author: ph-u
# script: eCology.sh
# desc: analysis ecology from simulation replicates
# in: bash eCology.sh [path/2/data]
# out: NA
# arg: 1 (req: data/*-{log,pri,sam}.csv)
# date: 20220509

[[ -z $1 ]] && grep -e "^\# desc\|^\# in" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit
pT=$1
for i in `ls ${pT}/*-log.csv`;do
	bNam=`echo -e ${i} | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "-"`
	tYpe=`echo -e ${bNam} | rev | cut -f 1 -d "_" | rev`
	rEp=`ls ${pT}/${bNam}*-sam.csv | wc -l`
	echo -e "analyzing: ${bNam} - (`date`)"
	Rscript eCology.r ${pT}/ ${bNam} ${tYpe} ${rEp}
done
exit
