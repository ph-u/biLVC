#!/bin/bash
# author: ph-u
# script: 00_master.sh
# desc: master pipeline for Approximate Bayesian Computing Rejection on time-series population dynamics
# in: bash 00_master.sh
# out: data/*-{par,pri,est,err,ext,nam}.txt ; data/fileList.txt ; raw/{customPrior,knownParam}.csv ; result/parDeduced.csv
# arg: 0
# date: 20211007,20211012,20211113 (fuse abcrLVC.sh),20211125

##### workbench prep #####
p0=`pwd`;p1=`dirname $0`
rEal=`echo -e "${p0}/${p1}" | sed -e "s/\/.$//"` #`realpath $0`
cd ${rEal} #`dirname ${rEal}`
sPhead="id,influencer,"
hEADer="${sPhead}value"
hEader="dateStamp,source,${hEADer},CI95Low,CI95High"
pRhead="${sPhead}par1,par2,type,min,max"
nUm='^[0-9]+$'
find .. | grep -e ".csv$" | grep -e "/raw/" | grep -e "_" > ../data/fileList.txt
[[ -f ../result/parDeduced.csv ]] || echo -e ${hEader} > ../result/parDeduced.csv
[[ -f ../raw/knownParam.csv ]] || echo -e ${hEader} > ../raw/knownParam.csv
[[ -f ../raw/customPrior.csv ]] || echo -e ${pRhead} > ../raw/customPrior.csv

##### pipeline capsule: one time-series one ABCR #####
while read -r lIne;do
	fIrst=`head -n 1 ${lIne}` ## header line listing all species names
	nSp=${fIrst//[^,]} ## count number of "," in header line
	nSp=$(( ${#nSp} +1 )) ## num of "," = num of id present (+1 for time column)
	nAm=`basename ${lIne} | cut -f 1 -d "."` ## name of time-series data file
	echo -e ${hEADer} > ../data/${nAm}-par.csv ## collect known parameters
	echo -e ${pRhead} > ../data/${nAm}-pri.csv ## collect unknown parameter priors

	## identify known & unknown parameters
	for i in `seq 2 ${nSp}`;do
		tArget=`echo -e ${fIrst} | cut -f ${i} -d ',' | tr -cd "[:print:]"`
		for j in r k `seq 2 ${nSp}`;do
			if [[ ${j} =~ ${nUm} ]];then # if ${j} is a number
				iNflu=`echo -e ${fIrst} | cut -f ${j} -d ',' | tr -cd "[:print:]"`
				mIn=-10
			else
				iNflu=${j}
				mIn=0
			fi
			if [[ ${j} == "k" ]];then mAx=25;else mAx=10;fi
			if [[ `grep -e ",${tArget},${iNflu}," ../raw/knownParam.csv | wc -l` -gt 0 ]];then
				grep -e ",${tArget},${iNflu}," ../raw/knownParam.csv | cut -f 3,4,5 -d "," >> ../data/${nAm}-par.csv ## collect self hindrance conditions
#			elif [[ `grep -e ",${tArget},${iNflu}," ../result/parDeduced.csv | wc -l` -gt 0 ]];then
#				grep -e ",${tArget},${iNflu}," ../result/parDeduced.csv | cut -f 3,4,5 -d "," >> ../data/${nAm}-par.csv ## collect known parameters from previous analyses
			elif [[ `grep -e "${tArget},${iNflu}," ../raw/customPrior.csv | wc -l` -gt 0  ]];then
				grep -e "${tArget},${iNflu}," ../raw/customPrior.csv | tail -n 1 >> ../data/${nAm}-pri.csv ## collect from custom-set priors; last entry as ref (in case of multi-entries)
			else
				echo -e "${tArget},${iNflu},0,0,Uniform,${mIn},${mAx}" >> ../data/${nAm}-pri.csv ## board define required prior
			fi;done;done

	## cluster run ABCR pipeline if needed
	if [[ $(( `wc -l < ../data/${nAm}-pri.csv` -1 )) -gt `grep -e ${nAm} ../result/parDeduced.csv | wc -l` ]];then
		echo -e "`date` >> ${nAm}"
		if [[ ${OSTYPE} == "linux-gnu" ]];then ## identify OS environment
			./slurm_submit.peta4-skylake_julia ${nAm} & # SLURM submission
			#julia abcrLVC.jl ${nAm} 1> ../data/${nAm}-est.txt 2> ../data/${nAm}-err.txt & # /dev/null
		else
			julia abcrLVC.jl ${nAm} 1> ../data/${nAm}-est.txt 2> ../data/${nAm}-err.txt # /dev/null
		fi
	else
		echo -e "`date` XX ${nAm}"
		rm ../data/${nAm}-pri.csv
	fi
done < ../data/fileList.txt

##### Wait until pipeline done #####
pP=9 # token to get detection loop started
while [[ ${pP} -gt 0 ]];do
	sleep 10
	if [[ ${OSTYPE} == "linux-gnu" ]];then ## identify OS environment
		pP=$[`gstatement -u pmh65 | grep -e "abcrLVC" | wc -l`-1]
	else
		pP=$[`ps aux | grep -e "abcrLVC" | wc -l`-1]
	fi
done

##### result restructure #####
cd ../data
ls *-est.txt > fileList-t.txt
while read -r lIne;do
	dAte=`date | cut -f 2,3,6 -d " "`
	nAm=`echo -e "${lIne}" | cut -f 1 -d "-"` # | rev | cut -f 1 -d "/" | rev`
	grep -e "^Parameter " ${nAm}-est.txt | cut -f 2 -d ":" | tr -d " \|)" | sed -e "s/(/,/g" > ${nAm}-ext.txt # cut -f 2 -d ":" | cut -f 1 -d "(" | tr -d " ""))"
	tail -n +2 ${nAm}-pri.csv | cut -f 1,2 -d "," | sed -e "s/^/${dAte},${nAm},/" > ${nAm}-nam.txt
	paste ${nAm}-nam.txt ${nAm}-ext.txt | sed -e "s/\t/,/" >> ../result/parDeduced.csv
	rm ${nAm}-{ext,nam}.txt
done < fileList-t.txt
#rm `ls | grep -e "\.txt$\|\.csv$" | grep -v "\-out"`
exit
