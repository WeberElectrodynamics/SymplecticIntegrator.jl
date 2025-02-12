
module WeberHamiltonian

using Symbolics: @variables, Differential, expand_derivatives, build_function
using Base.Iterators: partition
using LinearAlgebra: norm
using NaNMath

function kinetic_energy(ps, d)
	KE = 0
	for momenta in partition(ps, d)
		p2 = sum(momenta .^ 2)
		KE += p2 / 2
	end
	return KE
end

function coulomb_potential(qs, d, charges, N)
	@assert length(charges) == N

	qs_mat = reshape(qs, d, :)
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

function ampere_potential(qs, ps, d, charges, c)
	N = length(qs) ÷ d

	qs_mat = reshape(qs, d, :)
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

function weber_potential_energy(qs, ps, d, charges, N, c)
	return coulomb_potential(qs, d, charges, N) - ampere_potential(qs, ps, d, charges, c)
end

function weber_hamiltonian(t, qs, ps, params)
	@assert length(qs) == length(ps)

	d, charges, N, c = params

	KE = kinetic_energy(ps, d)
	PE = weber_potential_energy(qs, ps, d, charges, N, c)

	return KE + PE
end

function weber_vector_fields(H::Function, params)
	d, charges, N, c = params
	total_q = d * N

	@variables qs[1:total_q] ps[1:total_q]
	H_sym = H(0, qs, ps, params)

	diff_expr(expr, var) = expand_derivatives(Differential(var)(expr))

	qdot_syms = [diff_expr(H_sym, ps[i]) for i in eachindex(ps)]
	pdot_syms = [-diff_expr(H_sym, qs[i]) for i in eachindex(qs)]

	qdot_exprs = [build_function(qdot, qs, ps, expression = Val{true}) for qdot in qdot_syms]
	pdot_exprs = [build_function(pdot, qs, ps, expression = Val{true}) for pdot in pdot_syms]

	return qdot_exprs, pdot_exprs
end

function save_function_expressions(exprs::Vector{Expr}, filename::String)
	open(filename, "w") do io
		for expr in exprs
			# Write the expression as a string.
			println(io, string(expr))
			# Write a delimiter to separate expressions.
			println(io, "------")
		end
	end
end

function generate_and_save_weber_vector_fields(d, charges, N, c)
	params = (d, charges, N, c)
	qdot_exprs, pdot_exprs = weber_vector_fields(weber_hamiltonian, params)
	save_function_expressions(qdot_exprs, "weber_vector_fields/qdot_exprs.txt")
	save_function_expressions(pdot_exprs, "weber_vector_fields/pdot_exprs.txt")
end

function load_function_expressions(filename::String)
	exprs = Expr[]
	file_content = read(filename, String)
	# Split on the delimiter. This returns a list of strings.
	for part in split(file_content, "------")
		part = strip(part)
		if !isempty(part)
			push!(exprs, Meta.parse(part))
		end
	end
	return exprs
end

function load_weber_vector_fields()
	qdot_exprs = load_function_expressions("weber_vector_fields/qdot_exprs.txt")
	pdot_exprs = load_function_expressions("weber_vector_fields/pdot_exprs.txt")

	qdot_funcs = [eval(expr) for expr in qdot_exprs]
	pdot_funcs = [eval(expr) for expr in pdot_exprs]

	return qdot_funcs, pdot_funcs
end

end
