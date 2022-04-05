#!/bin/bash
#! author: ph-u
#! script: pipe.sh
#! desc: Bayesian Inference pipeline
#! in: bash pipe.sh
#! out: (SEE bayesInfer.r)
#! arg: 0
#! date: 20220101, 20220108 (main pipeline), 20220331 (+generalized Lotka-Volterra option)

##### set enviroment ##### 2022{0101,0103}
nUm1='^[0-9]';nUm2='[0-9]$' #nUm='^[0-9]+*[0-9]+$' # line start & end with numbers
p1=`pwd`;p2=`dirname $0`;p0=`echo -e "${p1}/${p2}"|sed -e "s/\.$//"`
echo -e "BImcmc-LV - (`date`)"

##### get time-series ##### 20220101
cd ${p0}
ls ../raw/*.csv | grep -v "p_" | rev | cut -f 2 -d "." | cut -f 1 -d "/" | rev > ../data/tmp # all csv basename
[[ -f ../data/fList.txt  ]] && rm ../data/fList.txt
while read -r L;do # get time-series data
	p=`head -n 2 ../raw/${L}.csv | tail -n 1`
	if [[ ${p} =~ ${nUm1} ]] && [[ ${p} =~ ${nUm2} ]];then # [[ ${p} =~ ${nUm} ]]
		[[ -f ../data/${L}-rec.txt ]] || echo -e ${L} >> ../data/fList.txt # list unprocessed time-series
	fi
done < ../data/tmp
rm ../data/tmp

##### Data Analysis ##### 2022{0101,0104,0108,0331}
while read -r L;do
	if [[ `echo ${L} | rev | cut -f 1 -d "_" | rev` == "LVC" ]];then
		tP="c";tP0="LVC"
	else
		tP="g";tP0="gLV"
	fi
	echo -e "${L} analysed by ${tP0} (`date`)"
	if [[ ${OSTYPE} == "linux-gnu" ]];then
		sbatch bayesInfer.r ${L} ${tP} 1> ../data/${L}-rec.txt &
	else
		Rscript bayesInfer.r ${L} ${tP} 1> ../data/${L}-rec.txt
	fi
done < ../data/fList.txt
rm ../data/fList.txt
echo -e "All queued - (`date`)"
##### sort program message after hpc done ##### 20220401
# for i in `ls ../data/*-rec.txt`;do j=`head -n 1 ${i} | rev | cut -f 1 -d " " | rev`;j0=`echo ${i} | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "-"`;cat slurm-${j}.out >> ../data/${j0}-rec.txt;done
exit
