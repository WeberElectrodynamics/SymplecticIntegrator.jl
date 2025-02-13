position(ϕ) = [cos(ϕ), sin(ϕ)]
λ = √(0.707)
velocity(ϕ) = [-λ * sin(ϕ), λ * cos(ϕ)]

function integrate(start, stop, step)
	phases = range(0, 2π, length = N + 1)

	particle1qs = position(phases[1])
	particle2qs = position(phases[2])
	particle3qs = position(phases[3])
	particle4qs = position(phases[4])

	particle1ps = velocity(phases[1])
	particle2ps = velocity(phases[2])
	particle3ps = velocity(phases[3])
	particle4ps = velocity(phases[4])

	q0 = [particle1qs..., particle2qs..., particle3qs..., particle4qs...]
	p0 = [particle1ps..., particle2ps..., particle3ps..., particle4ps...]

	tspan = (start, stop)
	tstep = step

	problem = GeometricIntegrators.HODEProblem(v!, f!, H, tspan, tstep, q0, p0)

	# integrator = GeometricIntegrator(problem, QinZhang())
	integrator = GeometricIntegrator(problem, ImplicitMidpoint())
	# integrator = GeometricIntegrator(problem, SRK3())
	# integrator = GeometricIntegrator(problem, Gauss(4))
	# integrator = GeometricIntegrator(problem, LobattoIIID(2))
	# integrator = GeometricIntegrator(problem, LobattoIIIE(2))
	# integrator = GeometricIntegrator(problem, LobattoIIIG(2))
	# integrator = GeometricIntegrator(problem, SymplecticEulerA())
	# integrator = GeometricIntegrator(problem, SymplecticEulerB())
	# integrator = GeometricIntegrator(problem, PartitionedGauss(1))

	sol = GeometricIntegrators.integrate(integrator)

	T = sol.t
	qs = sol.s.q
	ps = sol.s.p

	return T, qs, ps
end

function _plot(t, qs, ps)
	trajectories = [reshape(q, d, N)' for q in qs]
	momenta = [reshape(p, d, N)' for p in ps]

	plt1 = plot(title = "Particle Trajectories", xlabel = "x", ylabel = "y", aspect_ratio = :equal)
	for i in 1:N
		xs = [traj[i, 1] for traj in trajectories]
		ys = [traj[i, 2] for traj in trajectories]

		pxs = [momenta[i, 1] for momenta in momenta]
		pys = [momenta[i, 2] for momenta in momenta]

		scatter!(plt1, [xs[1]], [ys[1]], label = "", marker = :circle)
		quiver!(plt1, [xs[1]], [ys[1]], quiver = ([pxs[1]], [pys[1]]), label = "", color = :red)

		plot!(plt1, xs, ys, label = charges[i])
	end

	total_energy = collect([H(0, q, p, params) for (q, p) in zip(qs, ps)])
	println(describe(total_energy))
	# println("\n\nVariance: ", var(total_energy), "\n\n")


	KE = [WeberHamiltonian.kinetic_energy(p, d) for p in ps]
	PE1 = [WeberHamiltonian.coulomb_potential(q, d, charges, N) for q in qs]
	PE2 = [WeberHamiltonian.ampere_potential(q, p, d, charges, c) for (q, p) in zip(qs, ps)]

	plt2 = plot(t, KE, label = "Kinetic Energy", xlabel = "Time", ylabel = "Energy")
	plot!(t, PE1, label = "Coulomb Potential Energy")
	plot!(t, PE2, label = "Ampere Potential Energy")

	plt3 = plot(total_energy, label = "Total Energy", xlabel = "Time", ylabel = "Energy")

	p = plot(
		plt1,
		plt2,
		plt3,
		layout = (3, 1),
		size = (800, 1600),
	)

	display(p)
end

@time t, qs, ps = integrate(0.0, 100, 0.001)
@time _plot(t, qs, ps)
