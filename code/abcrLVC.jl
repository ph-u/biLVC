#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: abcrLVC.jl
# desc: ABCR-LVC Bayesian Inference on natural-logged time-series
# in: julia --threads 3 abcrLVC.jl [basename of time-series]
# out: [basename of time-series]-est.txt
# arg: 1
# date: 20211101, 20211121 (fuse logData.r)

using ApproxBayes, CSV
include("abcrSIM.jl")
nAm = ARGS[1]
pAth = ifelse(Sys.islinux(), "/home/pmh65/pj/0_abcrPar/", "../")

rawTS = CSV.read(pAth*"raw/"*nAm*".csv", DataFrame; header=1, type=Float64)
kNown = CSV.read(pAth*"data/"*nAm*"-par.csv", DataFrame; header=1, types=[String, String, Float64])
pRior = CSV.read(pAth*"data/"*nAm*"-pri.csv", DataFrame; header=1, types=[String, String, Float64, Float64, String, Float64, Float64])

##### population data log-transformation #####
for i=1:nrow(rawTS), j=2:ncol(rawTS)
	rawTS[i,j] = log(rawTS[i,j]+1)
end

CSV.write(pAth*"data/"*nAm*"-log.csv", rawTS)
prIor, cSt = exPar(pRior, kNown, rawTS)

##### ABCR-LVC #####
abcR = ABCRejection(dist4abc, nrow(prIor), .2, Prior([cDistrib(prIor[i,:]) for i in 1:nrow(prIor)]); maxiterations = Int(1e7), constants = cSt)

q = runabc(abcR, rawTS, parallel=true)
string(q) # print decision summary to stdout for bash collection
writeoutput(q; dir=pAth*"data/", file=nAm*"-out.txt") # output successful parameter combinations (default 100 in ABCRejection)
