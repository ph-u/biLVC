#!/bin/bash
# author: ph-u
# script: biLVary.sh
# desc: SLURM-optimized Bayesian Inference Lotka-Volterra multi-replicate pipeline
# in: bash biLVary.sh [prep (1) / BI-MCMC (2) / post (3)] [replicates] [max-iteration (in 10^5)] [raw data directory name (must be same level as src/)]
# out: NA
# arg: 4
# date: 20231109 (supersede biLVrep.sh & eCology.sh)

[[ -z $1 ]] && grep -e "^\# desc\|^\# in" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit
mkdir -p ../{data,result}
sT=$1
[[ -z $2 ]] && rP=1 || rP=$2
[[ -z $3 ]] && mX=1 || mX=$3
[[ -z $4 ]] && ptIN="raw" || ptIN=$4
p1=`pwd`;p2=`dirname $0` # get full path of pipeline
p3=`echo -e "${p1}/${p2}"|sed -e "s/\.$//"`
[[ ${sT} -gt 4 ]] && echo -e "Check input -- might be wrong?" && exit

echo -e "BImcmc-LV - (`date`)"
##### List parameters #####
[[ -f "../parIN.csv" ]]&& rm ../parIN.csv
for rAw in `ls ../${ptIN}/*.csv | grep -v "^p_" | rev | cut -f 2 -d "." | cut -f 1 -d "/" | rev`;do #time-series basename
    [[ `echo ${rAw} | rev | cut -f 1 -d "_" | rev` == "LVC" ]]&&tP="c"||tP="g"
    if [[ ${sT} -eq 1 ]];then
        echo -e "${rAw},${tP},${rP},${p1}" >> ../parIN.csv
    elif [[ ${sT} -eq 2 ]];then
        for rEp in `seq 1 ${rP}`;do echo -e "${rAw},${tP},${mX},${rEp},${p1}" >> ../parIN.csv;done
    else
        echo -e "${p3}/../data/,${rAw},`echo ${rAw} | rev | cut -f 1 -d "_" | rev`,${rP}" >> ../parIN.csv
    fi
done
echo -e "Parameter listing done - (`date`)"

##### Assemble run env #####
if [[ ${sT} -eq 1 ]];then
    nAm="tsFIT"
    sRt="tsProcess.r"
    tIme="12:00:00"
elif [[ ${sT} -eq 2 ]];then
    nAm="biLV"
    sRt="bayesInfer.r"
    tIme="12:00:00"
else
    nAm="ecolAna"
    sRt="eCology.r"
    tIme="00:10:00"
fi

##### Assemble run script #####
cat arrayChild.sh > aRrayC.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${nAm}\n#SBATCH --time=${tIme}\n#SBATCH --array=1-`wc -l < ../parIN.csv`\n\ncd ${p3}\nRscript ${sRt} ${p0}" >> aRrayC.sh
echo -e "Run script assemble done - (`date`)"

##### Execute #####
sbatch aRrayC.sh &
exit
