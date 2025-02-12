using GLMakie

fig = Figure()
ax = Axis(fig[1, 1])
x = range(0, 10π, length = 100)
phase = Observable(0.0)
y = @lift(sin.(x .+ $phase))
lines!(ax, x, y)
display(fig)

@async while isopen(fig.scene)
	phase[] += 0.1
	sleep(0.01)
end
