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
c = 10

function kinetic_energy(ps)
	KE = 0
	for momenta in partition(ps, d)
		p2 = sum(momenta .^ 2)
		KE += p2 / 2
	end
	return KE
end

function coulomb_potential(qs, charges)
	qs_mat = reshape(qs, d, N)
	qs_mat = transpose(qs_mat)

	PE = 0
	for i in 1:N
		for j in i+1:N
			rel_qs = qs_mat[i, :] - qs_mat[j, :]
			r = norm(rel_qs)
			PE += charges[i] * charges[j] / r
		end
	end

	return PE
end

function ampere_potential(qs, ps)
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

			PE += (charges[i] * charges[j] / r) * ṙ^2 / 2c^2
		end
	end

	return PE
end

function potential_energy(qs, ps, charges)
	return coulomb_potential(qs, charges) - ampere_potential(qs, ps)
end

function H(t, qs, ps, params)
	@assert length(qs) == d * N && length(ps) == d * N

	KE = kinetic_energy(ps)
	PE = potential_energy(qs, ps, charges)

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
p0 = [particle1ps..., particle2ps..., particle3ps..., particle4ps...]

tspan = (0.0, 21)
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
	plot!(plt1, xs, ys, label = charges[i])
	scatter!(plt1, [xs[1]], [ys[1]], label = "", marker = :circle)
end

total_energy = collect([H(t_i, q, p, nothing) for (t_i, q, p) in zip(t, qs, ps)])
using StatsBase
println(describe(total_energy))

KE = [kinetic_energy(p) for p in ps]
PE1 = [coulomb_potential(q, charges) for q in qs]
PE2 = [ampere_potential(q, p) for (q, p) in zip(qs, ps)]

plt2 = plot(t, KE, label = "Kinetic Energy", xlabel = "Time", ylabel = "Energy")
plot!(t, PE1, label = "Coulomb Potential Energy")
plot!(t, PE2, label = "Ampere Potential Energy")

plot(plt1, plt2, layout = (3, 1), size = (800, 1600), aspect_ratio = :equal)

