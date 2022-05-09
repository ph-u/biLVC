#!/bin/bash
# author: ph-u
# script: sortMsg.sh
# desc: sort CSD3 slurm output captures into record files
# in: bash sortMsg.sh
# out: (output captures sorting)
# arg: 0
# date: 20220508 (segregate from biLVrep.sh)

##### sort program message after hpc done ##### 20220401
for i in `ls ../data/*-rec.txt`;do
	j=`head -n 1 ${i} | rev | cut -f 1 -d " " | rev`
	cat slurm-${j}.out >> ../data/${i}
done
exit

