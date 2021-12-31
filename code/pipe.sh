#!/bin/bash
#! author: ph-u
#! script: 0_pipe.sh
#! desc: ABCnLVC pipeline
#! in: bash 0_pipe.sh
#! out: (SEE respective scripts)
#! arg: 0
#! date: 20211226

nAm=$1
sCript="0_analyse.jl"
if [[ ${OSTYPE} == "linux-gnu"  ]];then
	srun julia ${sCript} ${nAm}
else
	julia ${sCript} ${nAm}
fi
