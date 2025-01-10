using Test
using SymplecticIntegrator

@testset "Hamiltonian Validation Tests" begin
	# function to test
	validate_hamiltonian = SymplecticIntegrator._validate_hamiltonian

	# valid hamiltonian
	valid_hamiltonian(q::Vector{Float64}, p::Vector{Float64})::Float64 = 0.0
	@test validate_hamiltonian(valid_hamiltonian) == true

	# invalid hamiltonian
	invalid_h1() = 0.0
	@test_throws ArgumentError validate_hamiltonian(invalid_h1)

	invalid_h2(q::Float64) = 0.0
	@test_throws ArgumentError validate_hamiltonian(invalid_h2)

	invalid_h3(q::Vector{Float64}, p::Vector{Float64})::Int = 0
	@test_throws ArgumentError validate_hamiltonian(invalid_h3)
end
