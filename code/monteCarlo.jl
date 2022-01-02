#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: monteCarlo.jl
# desc: monte carlo functions
# in: include("monteCarlo.jl")
# out: none
# arg: 0
# date: 20211228
##reading cites
## https://doi.org/10.1111/j.2041-210X.2011.00179.x
## https://doi.org/10.1002/2014WR015386

include("LotkaVolterra.jl")
include("dataHandle.jl")

##### simulation monomer ##### 20211230
function runModel(pRi::DataFrame, dAta::DataFrame)
	p = [rand(cDistrib(pRi[i,:])) for i in 1:nrow(pRi)]
	x0 = [median(dAta[findall(dAta[:,1] .== min(dAta[:,1]...)),i]) for i in 2:ncol(dAta)]
	tT = sort(unique(dAta[:,1]))
	iNt = min([tT[i]-tT[i-1] for i in 2:length(tT)]...)
	z = ode!(x0,p,extrema(dAta[:,1]),iNt)
	z0 = wHich(z,1,dAta[:,1])
	return tsMatch(z0,dAta), p
end

##### Rejection algorithm ##### 20211230
## data structure: successful trial [row] * (parameter set, iteration, match percentage) [col]
function ABCr(dAta::DataFrame, pRi::DataFrame, aCc=.5, mAxIt=1e5, sAmple=100)
	rEc = DataFrame(zeros(0,nrow(pRi)+2), :auto)
	i, iTer = 0, 0 # tags
	while (i<sAmple) && (iTer<Int64(mAxIt))
		pC, pAra = runModel(pRi, dAta)
		if (pC>=aCc)
			i += 1
			push!(rEc, append!(pAra,[iTer,pC]))
		end
		iTer += 1
	end
	return rEc
end

##### Markov Chain Monte Carlo (MCMC) ##### 20211231
## MCMC: Metropolis–Hastings algorithm
function ABCmcmc(dAta::DataFrame, pRi::DataFrame, aCc=.95, mAxIt=1e5, sAmple=100, cHain=4)
	cHAins = ABCr(dAta,pRi,0.,mAxIt,sAmple) # chain start points
	cHains = wHich(cHAins,ncol(cHAins),reverse(sort(cHAins[:,ncol(cHAins)]))[1:cHain])
	pRiCP = deepcopy(pRi)
	rEc = DataFrame(zeros(0,nrow(pRi)+2), :auto)
	i, iTer, sTag, cH = 0, 0, 0, 1 # tags
	while (cH<=cHain)
		if (i == 0) ## modify prior
			pRiCP[:,3] = v(cHains[cH,1:nrow(pRi)])
			pRiCP[:,4] = (pRiCP[:,7] .- pRiCP[:,6]) ./ 100
			pRiCP[:,5] .= "Normal"
			pcLast = cHains[cH,ncol(cHains)]
		elseif (pC > pcLast)
			pRiCP[:,3] = pAra
			pcLast = min(pC,aCc)
		end
		pC, pAra = runModel(pRiCP, dAta)
		sTag = ifelse(pC-pcLast < .1, sTag+1, 0)
		if (pC>=aCc)
			i += 1
			pAcp = deepcopy(pAra)
			push!(rEc,append!(pAcp,[cH,pC]))
		end
		iTer += 1
		if (i>sAmple) || (iTer%Int64(mAxIt)==0) || (sTag>sAmple)
			i, sTag = 0, 0
			cH += 1
		end
	end
	return rEc
end
