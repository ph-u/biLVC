#!/bin/bash
#! author: ph-u
#! script: pipe.sh
#! desc: Bayesian Inference pipeline
#! in: bash pipe.sh
#! out: (SEE bayesInfer.r)
#! arg: 0
#! date: 20220101, 20220108 (main pipeline)

##### set enviroment ##### 2022{0101,0103}
nUm1='^[0-9]';nUm2='[0-9]$' #nUm='^[0-9]+*[0-9]+$' # line start & end with numbers
p1=`pwd`;p2=`dirname $0`;p0=`echo -e "${p1}/${p2}"|sed -e "s/\.$//"`
echo -e "BImcmc-LVC - (`date`)"

##### get time-series ##### 20220101
cd ${p0}
ls ../raw/*.csv | grep -v "p_" | rev | cut -f 2 -d "." | cut -f 1 -d "/" | rev > ../data/tmp # all csv basename
[[ -f ../data/fList.txt  ]] && rm ../data/fList.txt
while read -r L;do # get time-series data
	p=`head -n 2 ../raw/${L}.csv | tail -n 1`
	if [[ ${p} =~ ${nUm1} ]] && [[ ${p} =~ ${nUm2} ]];then # [[ ${p} =~ ${nUm} ]]
		[[ -f ../data/${L}-est.txt ]] || echo -e ${L} >> ../data/fList.txt # list unprocessed time-series
	fi
done < ../data/tmp
rm ../data/tmp

##### Data Analysis ##### 2022{0101,0104,0108}
while read -r L;do
	echo -e "${L} (`date`)"
	if [[ ${OSTYPE} == "linux-gnu" ]];then
		sbatch bayesInfer.r ${L} 1> ../data/${L}-est.txt &
	else
		Rscript bayesInfer.r ${L} 1> ../data/${L}-est.txt
	fi
done < ../data/fList.txt
#rm ../data/fList.txt
echo -e "All queued - (`date`)"
exit
