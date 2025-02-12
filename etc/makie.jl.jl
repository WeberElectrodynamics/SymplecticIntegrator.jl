using GLMakie

fig = Figure(backgroundcolor = :grey, size = (600, 600))

l = 1.2

ax1 = Axis(
	fig[1, 1],
	aspect = 1,
	xlabel = "x",
	ylabel = "y",
	limits = (-l, l, -l, l),
)

phase = Observable(0.0)
point = lift(p -> Point2(cos(p), sin(p)), phase)

scatter!(ax1, point)

fig

@async while true
	sleep(0.01)
	phase[] += 0.01
end
