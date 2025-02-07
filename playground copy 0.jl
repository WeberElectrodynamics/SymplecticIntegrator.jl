using Symbolics
using Base.Iterators: partition
using LinearAlgebra: norm
using GeometricIntegrators
using Plots

const d = 2
const N = 4
const total_q = d * N

@variables qs[1:total_q] ps[1:total_q]

charges = [1, -1, 1, -1]
c = 10e-3

function H(t, qs, ps, params)
	@assert length(qs) == d * N && length(ps) == d * N

	KE = 0
	for momenta in partition(ps, d)
		p2 = sum(momenta .^ 2)
		KE += p2 / 2
	end

	qs_mat = reshape(qs, d, N)
	qs_mat = transpose(qs_mat)

	momenta = collect(partition(ps, d))

	PE = 0
	for i in 1:N
		for j in i+1:N
			rel_qs = qs_mat[i, :] - qs_mat[j, :]

			ps_i = momenta[i]
			ps_j = momenta[j]

			rel_ps = ps_i - ps_j

			r = norm(rel_qs)
			ṙ = sum(rel_qs .* rel_ps) / r

			PE += charges[i] * charges[j] / r * (1 - ṙ^2 / 2c^2)
		end
	end

	return KE + PE
end

H_sym = H(0, qs, ps, 0)

diff_expr(expr, var) = expand_derivatives(Differential(var)(expr))

qdot_syms = [diff_expr(H_sym, ps[i]) for i in eachindex(ps)]
pdot_syms = [-diff_expr(H_sym, qs[i]) for i in eachindex(qs)]

qdot_funcs = [eval(build_function(qdot, qs, ps)) for qdot in qdot_syms]
pdot_funcs = [eval(build_function(pdot, qs, ps)) for pdot in pdot_syms]

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

particle1qs = [1.0, 0.0]
particle2qs = [0.0, 1.0]
particle3qs = [-1.0, 0.0]
particle4qs = [0.0, -1.0]

particle1ps = [0.0, 0.5]
particle2ps = [-0.5, 0.0]
particle3ps = [0.0, -0.5]
particle4ps = [0.5, 0.0]

q0 = [particle1qs..., particle2qs..., particle3qs..., particle4qs...]
# q0 = [particle1qs..., particle3qs...]
p0 = [particle1ps..., particle2ps..., particle3ps..., particle4ps...]
# p0 = [particle1ps..., particle3ps...]

tspan = (0.0, 5.0)
tstep = 0.001

prob = HODEProblem(v!, f!, H, tspan, tstep, q0, p0)

int = GeometricIntegrator(prob, Gauss(2))
sol = integrate(int)

t = sol.t
qs = sol.s.q
ps = sol.s.p

N_particles = length(q0) ÷ d
num_steps = length(t)

trajectories = [reshape(q, d, N_particles)' for q in qs]

plt1 = plot(title = "Particle Trajectories", xlabel = "x", ylabel = "y")
for i in 1:N_particles
	xs = [traj[i, 1] for traj in trajectories]
	ys = [traj[i, 2] for traj in trajectories]
	plot!(plt1, xs, ys, label = "Particle $i")
	scatter!(plt1, [xs[1]], [ys[1]], label = "", marker = :circle)
end

energies = [H(t_i, q, p, nothing) for (t_i, q, p) in zip(t, qs, ps)]

global_errors = [abs(E - energies[1]) for E in energies]
local_errors = [abs(energies[i+1] - energies[i]) for i in 1:(length(energies)-2)]

plt2 = plot(t, global_errors, title = "Global Energy Error", xlabel = "Time", ylabel = "Error", label = "Global")
plt3 = plot(t[2:end], local_errors, title = "Local Energy Error", xlabel = "Time", ylabel = "Error", label = "Local")

plot(plt1, plt2, plt3, layout = (3, 1), size = (400, 1200), aspect_ratio = :equal)
