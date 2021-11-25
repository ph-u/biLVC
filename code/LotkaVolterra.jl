#!/home/pmh65/bin/julia-1.6.3/bin/julia
# author: ph-u
# script: LotkaVolterra.jl
# desc: Lotka-Volterra scalable model
# in: jl-script `include("LotkaVolterra.jl")`
# out: none
# arg: 0
# date: 20210624, 20211104, 20211109

using DifferentialEquations

##### Lotka-Volterra ODE model #####
# log: number scale >> numeric value
# p: [r1,k1,a11->a1n,r2,k2,a21->a2n...]
function LV!(dx,x,p,t)
	x[x .< 0] .= 0
	nSp = length(x)
	for i in eachindex(x)
		w = nSp*(i-1)+1 ## start of appropriate parameter set
		q = x[i]*p[w] ## position of parameter in set
		dx[i] = deepcopy(q)
		for j in eachindex(x)
			dx[i] -= q*x[j]*p[w+1+j]/p[w+1]
		end
	end
	dx[x .== 0] .= 0
end

##### ODE model wrap #####
function ode!(x0::Vector{Float64}, p::Vector, t::Tuple{Float64, Float64}, iNterval::Float64)
	dE = ODEProblem(LV!, x0, t, p)
	return solve(dE, Feagin14(), saveat=iNterval) # sufficiently deterministic integration
end
