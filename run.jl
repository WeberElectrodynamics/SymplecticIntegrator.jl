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

tspan = (0.0, 10)
tstep = 0.001

problem = HODEProblem(v!, f!, H, tspan, tstep, q0, p0)
integrator = GeometricIntegrator(problem, Gauss(2))
sol = integrate(integrator)

t = sol.t
qs = sol.s.q
ps = sol.s.p

trajectories = [reshape(q, d, N)' for q in qs]
# momenta = [reshape(p, d, N_particles)' for p in ps]

plt1 = plot(title = "Particle Trajectories", xlabel = "x", ylabel = "y", aspect_ratio = :equal)
for i in 1:N
	xs = [traj[i, 1] for traj in trajectories]
	ys = [traj[i, 2] for traj in trajectories]
	plot!(plt1, xs, ys, label = charges[i])
	scatter!(plt1, [xs[1]], [ys[1]], label = "", marker = :circle)
end

total_energy = collect([H(0, q, p, params) for (q, p) in zip(qs, ps)])
println(describe(total_energy))

KE = [WeberHamiltonian.kinetic_energy(p, d) for p in ps]
PE1 = [WeberHamiltonian.coulomb_potential(q, d, charges, N) for q in qs]
PE2 = [WeberHamiltonian.ampere_potential(q, p, d, charges, c) for (q, p) in zip(qs, ps)]

plt2 = plot(t, KE, label = "Kinetic Energy", xlabel = "Time", ylabel = "Energy")
plot!(t, PE1, label = "Coulomb Potential Energy")
plot!(t, PE2, label = "Ampere Potential Energy")


p = plot(
	plt1,
	plt2,
	layout = (2, 1),
	size = (800, 1600),
)

display(p)
