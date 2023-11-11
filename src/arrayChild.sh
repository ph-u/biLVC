#!/bin/bash
# author: ph-u
# script: arrayChild.sh
# desc: get reassemble & call respective script to run
# in: bash arrayChild.sh [step 1/2/3] [slurm array index]
# out: NA
# arg: 2
# date: 20231110

#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p icelake-himem
