#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: LotkaVolterra.jl
# desc: Lotka-Volterra scalable model
# in: jl-script `include("LotkaVolterra.jl")`
# out: none
# arg: 0
# date: 20210624, 20211104, 20211109

using DifferentialEquations

##### parameter structure #####
struct param
	r::Vector{Float64}
	K::Vector{Float64}
	a::Array{Float64,2}
end

##### Lotka-Volterra ODE model #####
# log: number scale >> numeric value
function LV!(dx,x,p,t)
	x[x .< 0] .= 0
	for i in eachindex(x)
		q = x[i]*p.r[i]
		dx[i] = deepcopy(q)
		for j in eachindex(x)
			dx[i] -= q*x[j]*p.a[i,j]/p.K[i]
		end;end
	dx[x .== 0] .= 0
end

##### map vector into param structure #####
# [r1, k1, a11->a1n, r2, k2, a21->a2n ... rn, kn, an1->ann]
function paraMap(nSp::Int64, pAr::Vector{Float64})
	r, k, a = zeros(nSp), zeros(nSp), zeros(nSp,nSp)
	a0, a1 = 1,1 # row, col of "a"
	for i in eachindex(pAr)
		q = i%(nSp+2)
		if (q==1)
			r[a0] = deepcopy(pAr[i])
		elseif (q==2)
			k[a0] = deepcopy(pAr[i])
		else
			a[a0,a1] = deepcopy(pAr[i])
			if (q==0);a0+=1;a0=min(a0,nSp);a1=0;end
			a1+=1
		end;end
	return param(r,k,a)
end

##### ODE model wrap #####
function ode!(x0::Vector{Float64}, p::Vector, t::Tuple{Float64, Float64}, iNterval::Float64)
	dE = ODEProblem(LV!, x0, t, paraMap(length(x0),p))
	return DataFrame(solve(dE, Feagin14(), saveat=iNterval)) # sufficiently deterministic integration
end
