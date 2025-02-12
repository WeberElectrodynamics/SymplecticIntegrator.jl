using BenchmarkTools
using FixedPointNumbers
using LinearAlgebra

# Define fixed‑point types with various precisions:
const Fixed32  = Fixed{Int64, 32}    # 64-bit integer with 32 fractional bits
const Fixed64  = Fixed{Int128, 64}   # 128-bit integer with 64 fractional bits
const FixedBig = Fixed{BigInt, 128}   # Arbitrary-precision integer with 128 fractional bits

#--------------------------------------------------------------------------
# Function: benchmark_operations
#
# Performs a mix of arithmetic and transcendental operations.
#
# For Float64, we directly use native functions.
# For fixed‑point types, basic arithmetic works natively.
# Transcendental functions are computed via conversion to Float64.
#--------------------------------------------------------------------------
function benchmark_operations(T, n)
	if T == Float64
		a = T(1.2345)
		b = T(6.789)
		s = zero(T)
		for i in 1:n
			s += a + b         # addition
			s += a * b         # multiplication
			s += a / b         # division
			s += sin(a)        # sine
			s += cos(b)        # cosine
			s += sqrt(a)       # square root
			s += log(b)        # logarithm
			s += a * a         # squaring
		end
		return s
	else
		# For all fixed-point types:
		a = T(1.2345)
		b = T(6.789)
		s = zero(T)
		for i in 1:n
			s += a + b         # addition
			s += a * b         # multiplication
			s += a / b         # division
			# Transcendental functions are computed via conversion:
			s += T(sin(Float64(a)))
			s += T(cos(Float64(b)))
			s += T(sqrt(Float64(a)))
			s += T(log(Float64(b)))
			s += a * a         # squaring (multiplication)
		end
		return s
	end
end

#--------------------------------------------------------------------------
# Function: dot_benchmark
#
# Computes the dot product on vectors of length n.
#
# For Float64, the built‑in dot function is used.
# For fixed‑point types, we generate random numbers (converted to T)
# and compute the dot product in a manual loop.
#--------------------------------------------------------------------------
function dot_benchmark(T, n)
	if T == Float64
		a = rand(T, n)
		b = rand(T, n)
		return dot(a, b)
	else
		a = T.(rand(n))
		b = T.(rand(n))
		s = zero(T)
		@inbounds for i in 1:n
			s += a[i] * b[i]
		end
		return s
	end
end

# Number of iterations (or vector length)
n = 10^5

println("----------------------------------------------------")
println("Benchmarking Mixed Operations (loop with $n iterations)")
println("----------------------------------------------------")

println("\nFloat64:")
@btime result64 = benchmark_operations(Float64, $n)

println("\nFixedPoint (Fixed32, Fixed{Int64,32}):")
@btime resultFixed32 = benchmark_operations(Fixed32, $n)

println("\nFixedPoint (Fixed64, Fixed{Int128,64}):")
@btime resultFixed64 = benchmark_operations(Fixed64, $n)

# println("\nFixedPoint (FixedBig, Fixed{BigInt,128}):")
# @btime resultFixedBig = benchmark_operations(FixedBig, $n)

println("\n----------------------------------------------------")
println("Benchmarking Dot Product on vectors of length $n")
println("----------------------------------------------------")

println("\nFloat64:")
@btime dot_result64 = dot_benchmark(Float64, $n)

println("\nFixedPoint (Fixed32):")
@btime dot_resultFixed32 = dot_benchmark(Fixed32, $n)

println("\nFixedPoint (Fixed64):")
@btime dot_resultFixed64 = dot_benchmark(Fixed64, $n)

# println("\nFixedPoint (FixedBig):")
# @btime dot_resultFixedBig = dot_benchmark(FixedBig, $n)
