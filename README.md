# EcoAtom

[![DOI](https://zenodo.org/badge/450426957.svg)](https://zenodo.org/badge/latestdoi/450426957)

- [(https://doi.org/10.3389/fmicb.2023.1178131](https://doi.org/10.3389/fmicb.2023.1178131)
- [https://github.com/ph-u/biLVC](https://github.com/ph-u/biLVC)
- [https://github.com/ph-u/UKCFRegDataSorting](https://github.com/ph-u/UKCFRegDataSorting)
- [https://github.com/ph-u/UKCFRegDNase](https://github.com/ph-u/UKCFRegDNase)

## Overview

Analytics for pairwise ecological interactions of taxa categories (interactome) using time-series abundance records. Pipeline is built based on Bayesian Inference adaptive Markov Chain Monte Carlo with generalized Lotka-Volterra equation.

## Requirements

- R >= 3.3.3
- R pkg: FME(1.3.6.2), deSolve(1.30), rootSolve(1.8.2.3), coda(0.19.4)
- GNU bash (platform tested):
- 4.2.46(2)-release (x86\_64-redhat-linux-gnu)
- 5.1.16(1)-release (x86\_64-pc-linux-gnu)

## Updates

### EcoAtom

`2.0.0` - SLURM system optimization  
	> Backward compatible with ver 1 (BaMGS)  
	> array jobs submission

### BaMGS

`1.0.1` - Bugs fixed  
	> Filnames with "p_" string now segregated with prefix "p_"; prefix "p_" is for skipping respective data from processing)  
	> R 4.3.1 adaptation on dataframe handling (eCology.r)  
	> Fix scenario when no simulations matching data (eCology.r)  
	- Modified CPU cluster partition (skylake retired in late Aug-2023, new nodes added into icelake in the CSD3 cluster)  
`1.0.0` Initial release ([https://doi.org/10.3389/fmicb.2023.1178131](https://doi.org/10.3389/fmicb.2023.1178131))

## Usage

1. `Rscript cutTime.r [basename] [time...]`
	- optional step
	- indicate the overall time-series file in [basename]
	- indicate every timepoint in data need to be separated
	- generate series of files between consecutive indicated cut times
0. `bash biLVary.sh 1 [replicate] [sample (*100,000)] [raw data directory name]`
	- generate prior boundaries & starting values
	- generate simulation seed numbers
	- generate population time-series values (auto-identify whether raw data need log\_e transformation or not)
0. `bash biLVary.sh 2 [replicate] [sample (*100,000)] [raw data directory name]`
	- plot MCMC trace graphs
	- output best-fit 0.1% (min output = 100) parameter combinations as spreadsheets
0. `bash biLVary.sh 3 [replicate] [sample (*100,000)] [raw data directory name]`
	- group parameter combinations into pairwise ecological relationships
0. ```for i in `[path/2/\*-rec.txt]`;do echo -e "${i}"; bash SimDataPlot.sh ${i};done```
	- optional step
	- plot simulations on top of time-series data
	- calculate proportion of simulations matching data (+/- 10% tolerance for single data point time-series, 50% data range for data with replicates)
0. `Rscript ecoBars.r [path/2/data] [sampleFreq_basename] [replicate] [simulation_output_size] [log10?]`
	- optional step (R >= v4.0.0)
	- plot sample size of groups (log10 or raw frequency)
	- plot overview stacked bar plot of pairwise ecological relationship distributions across all replicates
