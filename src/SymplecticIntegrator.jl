module SymplecticIntegrator

using Symbolics

# Derive the Hamiltonian eqns and return them as Julia expressions
function _derive_hamiltonian_eqns(H::Function, dims)
	# Validate Hamiltonian
	_validate_hamiltonian(H)

	# Declare symbolic variables
	@variables qs[1:dims] ps[1:dims]

	# Declare symbolic expr
	H_ = H(qs, ps)

	# Declare differentiator
	diff(expr::Num, var::Num) = expand_derivatives(Differential(var)(expr))

	# Derive Hamiltonian eqns and build julia expressions
	q̇s = [
		eval(build_function(diff(H_, ps[i]), qs, ps))
		for i in 1:dims
	]

	ṗs = [
		eval(build_function(-diff(H_, qs[i]), qs, ps))
		for i in 1:dims
	]

	return [q̇s, ṗs]
end

# Validation function for the user provided Hamiltonian
function _validate_hamiltonian(H::Function)
	# Check H has only 1 method
	if length(methods(H)) != 1
		throw(ArgumentError("Hamiltonian must have only one method."))
	end

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
	if !(arg_types[1] <: AbstractVector{<:Number} && arg_types[2] <: AbstractVector{<:Number})
		throw(ArgumentError("Hamiltonian must accept two AbstractVector{<:Number} arguments."))
	end

	# Get return type (use concrete arg types so Julia compiler can infer return type)
	return_types = Base.return_types(H, (Vector{Float64}, Vector{Float64}))
	return_type = return_types[1]

	# Check return type
	if !(return_type <: Number)
		throw(ArgumentError("Hamiltonian must return a Number."))
	end

	return true
end

end
