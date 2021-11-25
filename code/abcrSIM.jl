#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: abcrSIM.jl
# desc: ABCR-LVC simulation + data matching main-function
# in: jl-script `include("abcrSIM.jl")`
# out: none
# arg: 0
# date: 20211101, 20211105, 20211108

using Statistics
include("LotkaVolterra.jl")
include("dataHandle.jl")

##### parameter values extraction #####
function exPar(pRior::DataFrame, kNown::DataFrame, rAw::DataFrame)
	nSp = ncol(rAw)-1
	nPar = (2+nSp)*nSp
	spNam = names(rAw)[2:end]
	iNflu = deepcopy(spNam)
	for pp in ["k","r"]; insert!(iNflu,1,pp);end
	pAr = DataFrame(id=repeat(spNam, inner=2+nSp), influ=repeat(iNflu, outer=nSp), prior=repeat([-999],nPar), known=repeat([-999.],nPar))

## fill either known / prior column
	for i in 1:nrow(pAr)
		prTMP = findall((pRior[:,1] .== pAr[i,"id"]) .& (pRior[:,2] .== pAr[i,"influ"]))
		knTMP = findall((kNown[:,1] .== pAr[i,"id"]) .& (kNown[:,2] .== pAr[i,"influ"]))
		if (length(knTMP) == 0)
			pAr[i,3] = prTMP[length(prTMP)]
		else
			pAr[i,4] = median(skipmissing(kNown[knTMP,3]))
		end
	end

	oUt = pAr[findall(pAr[:,"prior"] .!= -999),"prior"] # vector of row numbers sequence in pRior
	return pRior[oUt,:], pAr[:,"known"] ## prior dataframe, constant vectors
end

##### matching data distance #####
function dist4abc(pRi::Vector, cSt::Vector, rawData)
	z = abcrMAP(pRi, cSt, rawData) #with known parameter values
	matchRatio(z, rawData), 1
end

##### standardize timesteps in simulation with data #####
function abcrMAP(pR0::Vector{Float64}, kN0::Vector, rAw::DataFrame)
	tm = sort(unique(rAw[:,1])) # timestep in raw data
	tP = rAw[(rAw[:,1] .== tm[1]), 2:end] # get initial conditions for all experimental replicates
	x0 = [median(skipmissing(tP[:,i])) for i in 1:ncol(tP)] # get simulation initial conditions
	s0 = min([tm[i] - tm[i-1] for i in 2:length(tm)]...) # min stepsize in data for simulation
	p0 = deepcopy(kN0)
	p0[findall(p0 .== -999.)] = pR0 ## fill in numbers from prior to known vector
	z = ode!(x0, kN0, extrema(tm), s0) # simulation

## data-matching
	if extrema(z.t) == extrema(tm) # ode truly successful
		z = DataFrame(z)
		sIm = z[[findall(x -> in(x,i), z[:,1])[1] for i in tm],:] # match sim-data time with raw
	else
		sIm = DataFrame(zeros(length(tm), length(x0)+1), Symbol.(names(rAw)))
	end
	return sIm
end
