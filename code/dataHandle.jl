#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: dataHandle.jl
# desc: data-reorganising functions
# in: jl-script `include("dataHandle.jl")`
# out: none
# arg: 0
# date: 20210624

using DataFrames, Distributions, StatsBase, Statistics

##### random numbers from indicated distribution ##### 20210624, 20210809
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

##### simulated value match ##### 20211228
## linear match ratio drop away from data range
function simMatch(sIm, dAta::Vector)
	x = extrema(dAta)
	if (x[1]<=sIm<=x[2])
		return 1
	elseif ((x[1]-1)<=sIm<=(x[2]+1))
		d = ifelse(sIm<x[1],x[1],x[2])
		return abs(sIm-d)
	else
		return -1
#		return max(1-abs(sIm/d-1),0) ## match range = [0,1]
	end
end

##### time-series match ##### 20211229
function tsMatch(tS::DataFrame, dAta::DataFrame)
	tIme = sort(unique(dAta[:,1])) ## data timeline
	if (nrow(tS) < length(tIme));return 0;end
	mAtch = 0 ## match ratio collection
	for i0 in 2:length(tIme)
		q1 = v(tS[findall(tS[:,1] .== tIme[i0]),2:ncol(tS)])
		q2 = dAta[findall(dAta[:,1] .== tIme[i0]),2:ncol(dAta)]
		for i1 in eachindex(q1)
			mAtch += simMatch(q1[i1],q2[:,i1])
		end;end
	return mAtch/((length(tIme)-1)*(ncol(tS)-1))
end

##### extract dataframe using reference column by vector ##### 20211230
function wHich(df::DataFrame, colNum::Int64, rEf::Vector)
	q = []
	for i in eachindex(rEf)
		append!(q,findall(df[:,colNum] .== rEf[i]))
	end
	return df[sort(unique(q)),:]
end

##### vectorize DataFrame col/row ##### 20211230
v = z -> vec(Array(z))
