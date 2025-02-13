using GeometricIntegrators
using Plots
using StatsBase

include("lib/weber_hamiltonian.jl")

d = 2
charges = [1, -1, 1, -1]
N = 4
c = 10

params = (d, charges, N, c)

# WeberHamiltonian.generate_and_save_weber_vector_fields(d, charges, N, c)

qdot_funcs, pdot_funcs = WeberHamiltonian.load_weber_vector_fields()
H = WeberHamiltonian.weber_hamiltonian

function v!(dq, t, q, p, params)
	for i in eachindex(dq)
		dq[i] = qdot_funcs[i](q, p)
	end
end

function f!(dp, t, q, p, params)
	for i in eachindex(dp)
		dp[i] = pdot_funcs[i](q, p)
	end
end
