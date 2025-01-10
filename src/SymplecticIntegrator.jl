module SymplecticIntegrator

# Validation function for the user provided Hamiltonian
function _validate_hamiltonian(H::Function)
	# Get first method
	H_ = methods(H)[1]

	# Get arg types
	arg_types = fieldtypes(H_.sig)

	# Check tuple length
	if length(arg_types) != 3
		throw(ArgumentError("Hamiltonian must accept two arguments."))
	end

	# Get arg types
	arg_types = arg_types[2:end] # First arg is typeof(H)

	# Check arg types
	if arg_types != (Vector{Float64}, Vector{Float64})
		throw(ArgumentError("Hamiltonian must accept two arrays of numeric types."))
	end

	# Get return type
	return_types = Base.return_types(H, (Vector{Float64}, Vector{Float64}))
	return_type = return_types[1]

	# Check return type
	if !(return_type == Float64)
		throw(ArgumentError("Hamiltonian must return a numeric Float64."))
	end

	return true
end

end
