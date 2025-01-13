using Test
using SymplecticIntegrator

# valid hamiltonians
function valid_hamiltonian(q::AbstractVector{<:Number}, p::AbstractVector{<:Number})::Number
	sum = 0
	for i in eachindex(q)
		sum += q[i] * p[i]
	end
	return sum
end


# invalid hamiltonians
function invalid_hamiltonian(q::AbstractVector{<:Number}, p::AbstractVector{<:Number}, x::Number)::Number
	sum = 0
	for i in eachindex(q)
		sum += q[i] * p[i] * x
	end
	return sum
end

@testset "Hamiltonian Validation Tests" begin
	# function to test
	validate_hamiltonian = SymplecticIntegrator._validate_hamiltonian

	# valid hamiltonians
	@test validate_hamiltonian(valid_hamiltonian) == true

	# invalid hamiltonians
	@test_throws ArgumentError validate_hamiltonian(invalid_hamiltonian)
end

@testset "Hamiltonian Eqns Derivation Tests" begin
	# function to test
	derive_hamiltonian_eqns = SymplecticIntegrator._derive_hamiltonian_eqns

	# dims
	dims = 3

	# valid hamiltonian
	result = derive_hamiltonian_eqns(valid_hamiltonian, dims)
	@test length(result[1]) == 3
	@test length(result[2]) == 3

	# invalid hamiltonian
	@test_throws ArgumentError derive_hamiltonian_eqns(invalid_hamiltonian, dims)
end
