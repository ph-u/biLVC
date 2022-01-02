#!/home/pmh65/bin/julia-1.6.3/bin/julia
#! author: ph-u
#! script: analyze.jl
#! desc: ABCnLVC main process
#! in: julia analyze.jl [time-series] [algorithm] [acceptance]
#! out: [time-series]{-[acceptance]}-[algorithm].csv
#! arg: 3
#! date: 20220101

#SBATCH -J abcnLVC-CPU
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07:00:00
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH -p skylake

##### environment #####
nAm,aLg,aCc = ARGS
pAth = ifelse(Sys.islinux(), "/home/pmh65/rds/0_abcrPar/", "../")

##### import #####
using CSV
include("monteCarlo.jl")
q = CSV.read(pAth*"data/"*nAm*"-log.csv", DataFrame; header=1, types=Float64)
p = CSV.read(pAth*"data/"*nAm*"-pri.csv", DataFrame; header=1, types=[String, String, Float64, Float64, String, Float64, Float64])

##### Approximate Bayesian Computing #####
if (aLg=="mcmc")
	q0 = ABCmcmc(q,p,parse(Float64,aCc))
else
	q0 = ABCr(q,p,parse(Float64,aCc))
end
CSV.write(pAth*"data/"*nAm*aCc*"-"*aLg*".csv", q0)

##### summary statistics #####
s = DataFrame(id=String[], influencer=String[], min=Float64[], low95CI=Float64[], Q1=Float64[], Q2=Float64[], Q3=Float64[], high95CI=Float64[], max=Float64[])
for i in 1:nrow(p)
	q1 = quantile(q0[:,i], [.05,.25,.5,.75,.95])
	q2 = extrema(q0[:,i])
	push!(s,[p[i,1],p[i,2],q2[1],q1[1],q1[2],q1[3],q1[4],q1[5],q2[2]])
end
CSV.write(pAth*"result/"*nAm*aCc*"-"*aLg*".csv", s)
