#!/home/pmh65/bin/julia-1.6.3/bin/julia
#! author: ph-u
#! script: 0_analyze.jl
#! desc: ABCnLVC main process
#! in: julia 0_analyze.jl [time-series basename]
#! out: [time-series basename]-{log,est}.csv
#! arg: 1
#! date: 20211226

#SBATCH -J abcnLVC-CPU
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07:00:00
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH -p skylake

using CSV
##include("abcrSIM.jl") #modify
aCc, nAm = ARGS
pAth = ifelse(Sys.islinux(), "/home/pmh65/rds/0_abcrPar/", "../")

rawTS = CSV.read(pAth*"raw/"*nAm*".csv", DataFrame; header=1, types=Float64)
#kNown = CSV.read(pAth*"data/"*nAm*"-par.csv", DataFrame; header=1, types=[String, String, Float64])
#pRior = CSV.read(pAth*"data/"*nAm*"-pri.csv", DataFrame; header=1, types=[String, String, Float64, Float64, String, Float64, Float64])

##### time-series log-transformation #####
for i=1:nrow(rawTS), j=2:ncol(rawTS)
	rawTS[i,j] = log(rawTS[i,j]+1)
end

CSV.write(pAth*"data/"*nAm*"-log.csv", rawTS)

##### initialize #####
a0 = ifelse(aCc=="r", .95, .5) # rej / mcmc
## MCMC: Metropolis–Hastings algorithm
