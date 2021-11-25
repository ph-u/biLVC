# Approximate Bayesian Computing Rejection on Lotka-Volterra Competition Model (ABCR-LVC) draft description

## Brief pipeline summary

A Bayesian Inference pipeline for time-series microbial dynamics. The master script `00\_master.sh` checks what time-series data needed to be analysed by the ABCR-LVC model; it collects known parameters for each time-series and pushes the necessary analyses in parallel into the [CSD3](https://www.hpc.cam.ac.uk/high-performance-computing) cluster.

## Important directories (d) / files (f)

name | type | notes
--- | --- | ---
`code/` | d | home for ABCR-LVC pipeline code files
`raw/` | d | manually-cleaned time-series files should be put here
`raw/customPrior.csv` | f | mannually-defined priors storage site
`result/parDeduced.csv` | f | date-stamped parameter values storage site

## Intermediate files tag (-[xxx]) table

tag | description
--- | ---
log | natural log of raw time-series data
par | known parameter collection for ABCR-LVC
pri | priors collection for ABCR-LVC
est | parameter estimation raw result
ext | temporary storage for estimated parameter values
nam | temporary storage for metadata for estimated parameter values

## Computing dependencies [language (L), package (P)]

Dependency | type | version
--- | --- | ---
bash | L | 4.2.46(2)-release
R | L | 3.3.3
julia | L | 1.6.3
ApproxBayes | P (julia) | 0.3.2
CSV | P (julia) | 0.9.5
DataFrames | P (julia) | 1.2.2
DifferentialEquations | P (julia) | 6.18.0
Distributions | P (julia) | 0.23.12
StatsBase | P (julia) | 0.33.10

## Conditions justifications

- self-hindrance values pre-set at 1: act as reference values for the equation subject; other interaction coefficients are compared relative to the intra-strain self-hindrance competition
