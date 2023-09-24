#!/bin/bash
# author: ph-u
# script: biLVrep.sh
# desc: multi-replicates Bayesian Inference Lotka-Volterra pipeline
# in: bash biLVrep.sh [prior prep (1) / BI-MCMC (2)] [replicates] [max-iteration (in 10^5)]
# out: data/*-rec.txt
# arg: 3
# date: 20220508 (supersede pipe.sh)

##### env #####
mkdir -p ../{data,result}
[[ -z $1 ]] && sT=1 || sT=$1
[[ -z $2 ]] && rP=1 || rP=$2
[[ -z $3 ]] && mX=1 || mX=$3
p1=`pwd`;p2=`dirname $0` # get full path of pipeline
echo -e "BImcmc-LV - (`date`)"
[[ ${sT} -gt 2 ]] && echo -e "Check input -- might be wrong?" && exit

##### analysis #####
cd `echo -e "${p1}/${p2}"|sed -e "s/\.$//"`
for rAw in `ls ../raw/*.csv | grep -v "^p_" | rev | cut -f 2 -d "." | cut -f 1 -d "/" | rev`;do # time-series basename
## basename id: LVC / gLV
	if [[ `echo ${rAw} | rev | cut -f 1 -d "_" | rev` == "LVC" ]];then
		tP="c";tP0="LVC"
	else
		tP="g";tP0="gLV"
	fi
## time-series process, priors calculation
	if [[ ${sT} == 1 ]];then
		if [[ ${OSTYPE} == "linux-gnu" ]];then
			sbatch tsProcess.r ${rAw} ${tP} ${rP} `pwd` 1> ../data/${rAw}-tsP.txt &
		else
			Rscript tsProcess.r ${rAw} ${tP} ${rP}
		fi
	fi
## Bayesian MCMC sampling with independent replicates
	if [[ ${sT} == 2 ]];then
		echo -e "${rAw}: eq ${tP0}, rep ${rP} (`date`)"
		for rEp in `seq 1 ${rP}`;do # replicates
			if [[ ${OSTYPE} == "linux-gnu" ]];then
				sbatch bayesInfer.r ${rAw} ${tP} ${mX} ${rEp} `pwd` 1> ../data/${rAw}-${rEp}-rec.txt &
			else
				Rscript bayesInfer.r ${rAw} ${tP} ${mX} ${rEp} 1> ../data/${rAw}-${rEp}-rec.txt
			fi
		done
	fi
done
echo -e "All queued - (`date`)"
exit
