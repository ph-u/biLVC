#!/bin/bash
#! author: ph-u
#! script: pipe.sh
#! desc: ABCnLVC pipeline
#! in: bash pipe.sh [algorithm] [acceptance]
#! out: (SEE respective scripts)
#! arg: 0
#! date: 20220101

##### set enviroment #####
nUm1='^[0-9]';nUm2='[0-9]$' #nUm='^[0-9]+*[0-9]+$' # line start & end with numbers
p1=`pwd`;p2=`dirname $0`;p0=`echo -e "${p1}/${p2}"|sed -e "s/\.$//"`
if [[ -z $1 ]] || [[ $1 == "mcmc" ]];then a="mcmc";else a="r";fi
[[ -z $2 ]] && t=.95 || t=$2
echo -e "ABCnLVC - ${a} algorithm at ${t} acceptance (`date`)"

##### get time-series #####
cd ${p0}
ls ../raw/*.csv | grep -v "p_" | rev | cut -f 2 -d "." | cut -f 1 -d "/" | rev > ../data/tmp # all csv basename
while read -r L;do # get time-series data
	p=`head -n 2 ../raw/${L}.csv | tail -n 1`
	if [[ ${p} =~ ${nUm1} ]] && [[ ${p} =~ ${nUm2} ]];then # [[ ${p} =~ ${nUm} ]]
		[[ -f ../data/${L}${t}-${a}.csv ]] || echo -e ${L} >> ../data/fList.txt # list unprocessed time-series
	fi
done < ../data/tmp
rm ../data/tmp

##### Approximate Bayesian Computing #####
while read -r L;do
	echo -e "${L} (`date`)"
	Rscript paraEst.r ${L} # parameter estimation
	if [[ ${OSTYPE} == "linux-gnu" ]];then
		srun julia analyze.jl ${L} ${a} ${t} 2> ${p0}../data/${L}-err.txt &
	else
		julia analyze.jl ${L} ${a} ${t} 2> ../data/${L}-err.txt
	fi
done < ../data/fList.txt
rm ../data/fList.txt
echo -e "ABCnLVC - ${a} algorithm at ${t} acceptance (main script done: `date`)"
exit
