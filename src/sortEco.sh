#!/bin/bash
# author: ph-u
# script: sortEco.sh
# desc: sort SWLURM output message + analysis ecology (cluster use)
# in: bash sortEco.sh [path/2/data]
# out: NA
# arg: 1 (req: data/*-{log,pri,sam}.csv, slurm-*.out)
# date: 20220510

[[ -z $1 ]] && grep -e "^\# desc\|^\# in" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit 

bash sortMsg.sh &
bash eCology.sh $1
exit
