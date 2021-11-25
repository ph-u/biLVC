#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: dataHandle.jl
# desc: data-reorganising functions
# in: jl-script `include("dataHandle.jl")`
# out: none
# arg: 0
# date: 20210624, 20210809, 20211101

using DataFrames, Distributions, StatsBase

##### random numbers from indicated distribution #####
function cDistrib(x)
	if x[5] == "Uniform" # low, up
		return Uniform(x[6],x[7])
	elseif in(x[5], ["Normal","Gaussian"]) # mu, sd
		return TruncatedNormal(x[3],x[4],x[6],x[7])
	elseif x[5] == "LogNormal" # mu, sd
		return truncated(LogNormal(x[3],x[4]),x[6],x[7])
	elseif x[5] == "Beta" # alpha, beta
		return truncated(Beta(x[3],x[4]),x[6],x[7])
	elseif x[5] == "Gamma" # alpha, theta
		return truncated(Gamma(x[3],x[4]),x[6],x[7])
	elseif x[5] == "Pareto" # alpha, theta
		return truncated(Pareto(x[3],x[4]),x[6],x[7])
	elseif x[5] == "InverseGamma" # alpha, theta
		return truncated(InverseGamma(x[3],x[4]),x[6],x[7])
	elseif x[5] == "GeneralizedExtremeValue" # mu, sigma ; shape = 0 (support IR)
		return truncated(GeneralizedExtremeValue(x[3],x[4],0),x[6],x[7])
	else
		return 1
	end
end

##### check proportion of simulation within time-stamped data range #####
function matchRatio(sIm::DataFrame, rAw::DataFrame)
	if sum(Array(sIm)) == 0;return 1;end # skip failed simulations
	pcMatch = 0
	for i in 1:nrow(sIm)
		t0 = rAw[findall(x -> in(x,sIm[i,1]), rAw[:,1]),:]
		for j in 2:ncol(t0)
			t1 = extrema(t0[:,j]) # quantile(t0[:,j], [.05,.95])
			pcMatch += t1[1] <= sIm[i,j] <= t1[2]
		end
	end
	return 1-pcMatch/(nrow(sIm)*(ncol(sIm)-1))
end
