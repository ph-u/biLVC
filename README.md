# BaMGS: *B*ayesian Inference *a*daptive *M*arkov Chain Monte Carlo using the *g*eneralized Lotka-Volterra equation and *s*tatistics framework

## Overview

Ecological analysis from time-series poly-category population size data to pairwise ecological reltionship distributions graphs.  Statistics can be performed on intermediate files.

## Framework stepwise commands

1. `Rscript cutTime.r [basename] [time...]`
	- optional step
	- indicate the overall time-series file in [basename]
	- indicate every timepoint in data need to be separated
	- generate series of files between consecutive indicated cut times
0. `bash biLVrep.sh 1 [replicate] [sample (*10,000)]`
	- generate prior boundaries & starting values
	- generate simulation seed numbers
	- generate population time-series values (auto-identify whether raw data need log\_e transformation or not)
0. `bash biLVrep.sh 2 [replicate] [sample (*10,000)]`
	- plot MCMC trace graphs
	- output best-fit 1% (min output = 100) parameter combinations as spreadsheets
0. `bash sortEco.sh [path/2/data] >> ../sortEco.log`
	- pair cluster out-message with cluster job number
	- group parameter combinations into pairwise ecological relationships
0. ```for i in `[path/2/\*-rec.txt]`;do echo -e "${i}"; bash SimDataPlot.sh ${i};done```
	- optional step
	- plot simulations on top of time-series data
	- calculate proportion of simulations matching data (+/- 10% tolerance for single data point time-series, 50% data range for data with replicates)
0. `Rscript ecoBars.r [path/2/data] [sampleFreq_basename] [replicate] [simulation_output_size] [log10?]`
	- optional step (R >= v4.0.0)
	- plot sample size of groups (log10 or raw frequency)
	- plot overview stacked bar plot of pairwise ecological relationship distributions across all replicates
