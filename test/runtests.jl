using CamparyBindings
using Test
using Random

const C = CamparyBindings

setprecision(BigFloat, 4096)

exact(x::C.CamparyFloat) = sum(BigFloat.(x.val))
prec(::Type{Float64}) = 53
prec(::Type{Float32}) = 24
# The theorems assume no underflow; a Float32 expansion of more than four limbs of
# p = 24 bits reaches the subnormal range.
maxlimbs(::Type{Float64}) = 8
maxlimbs(::Type{Float32}) = 4

# A random ulp-nonoverlapping expansion with K limbs, |x_{i+1}| <= ulp(x_i), with leading
# limb of size about 2^e.
function rand_expansion(::Type{T}, K; e = 0) where {T}
	v = zeros(T, K)
	v[1] = T(ldexp(rand() + 0.5, e)) * rand((-1, 1))
	for i in 2:K
		v[i] = rand(T) * eps(v[i-1]) * rand((-1, 1))
	end
	return C.CamparyFloat{K}(v)
end

const is_ulp_nonoverlapping = C.is_ulp_nonoverlapping

Random.seed!(20261126)

@testset "CamparyBindings.jl" begin

@testset "error-free transforms ($T)" for T in (Float64, Float32)
	for _ in 1:1000
		a, b = T(randn() * 2.0^rand(-30:30)), T(randn() * 2.0^rand(-30:30))
		s, e = TwoSum(a, b)
		@test BigFloat(s) + BigFloat(e) == BigFloat(a) + BigFloat(b)
		hi, lo = abs(a) >= abs(b) ? (a, b) : (b, a)
		s, e = FastTwoSum(hi, lo)
		@test BigFloat(s) + BigFloat(e) == BigFloat(hi) + BigFloat(lo)
		s, e = FastTwoDiff(hi, lo)
		@test BigFloat(s) + BigFloat(e) == BigFloat(hi) - BigFloat(lo)
		s, e = C.TwoDiff(a, b)
		@test BigFloat(s) + BigFloat(e) == BigFloat(a) - BigFloat(b)
		p, e = TwoProd(a, b)
		@test BigFloat(p) + BigFloat(e) == BigFloat(a) * BigFloat(b)
	end
end

@testset "VecSum and FastVecSum are error-free ($T)" for T in (Float64, Float32)
	# regression: for K = 2 the first limb used to be added twice
	x = C.CamparyFloat{2}(T[1, ldexp(T(1), -40)])
	@test exact(VecSum(x)) == exact(x)
	@test VecSum(C.CamparyFloat{1}(T[3])).val == T[3]
	for K in 1:maxlimbs(T), _ in 1:200
		# VecSum is error-free on arbitrary input
		x = C.CamparyFloat{K}(T.(randn(K) .* 2.0 .^ rand(-40:40, K)))
		@test exact(VecSum(x)) == exact(x)
		r = C.CamparyFloat{K}(zeros(T, K))
		@test exact(VecSum(x, r)) == exact(x)
		y = C.CamparyFloat{K}(copy(x.val))
		@test exact(VecSum!(y)) == exact(x)
		@test y.val == VecSum(x).val
		# FastVecSum needs the exponent condition of Fast2Sum: decreasing input
		x = rand_expansion(T, K)
		@test exact(FastVecSum(x)) == exact(x)
		y = C.CamparyFloat{K}(copy(x.val))
		@test exact(FastVecSum!(y)) == exact(x)
		@test FastVecSum(x).val == VecSum(x).val
	end
end

@testset "VecSum4Add1 ($T)" for T in (Float64, Float32)
	for K in 1:maxlimbs(T)-1, _ in 1:100
		x = rand_expansion(T, K); y = T(randn())
		r = C.VecSum4Add1(x, y)
		@test length(r) == K + 1
		@test exact(r) == exact(x) + BigFloat(y)
		@test C.VecSum4Add1(x, y, C.CamparyFloat{K+1}(zeros(T, K+1))).val == r.val
	end
	@test_throws DimensionMismatch C.VecSum4Add1(rand_expansion(T, 3), one(T), C.CamparyFloat{3}(zeros(T, 3)))
end

@testset "VecSumErrBranch and Renorm ($T)" for T in (Float64, Float32)
	# regression: the early return used to hand back a Vector
	x = C.CamparyFloat{3}(T[1, ldexp(T(1), -2prec(T)), ldexp(T(1), -4prec(T))])
	f = C.VecSumErrBranch(x, Val(1))
	@test f isa C.CamparyFloat{T, 1}
	@test f.val == T[1]
	# regression: the in-place form used to keep stale limbs in its output
	r = C.CamparyFloat{4}(fill(T(7), 4))
	@test C.VecSumErrBranch(C.CamparyFloat{2}(T[1, 0]), r).val == T[1, 0, 0, 0]
	@test C.FastVecSumErrBranch(C.CamparyFloat{2}(T[1, 0]), C.CamparyFloat{4}(fill(T(7), 4))).val == T[1, 0, 0, 0]

	for K in 2:maxlimbs(T)+2, _ in 1:200
		# input that overlaps by at most p-2 digits, in decreasing order: the merge of two
		# ulp-nonoverlapping expansions, as in Theorem 3.3.2 of Popescu's thesis
		a, b = rand_expansion(T, K ÷ 2 + 1), rand_expansion(T, K - K ÷ 2; e = rand(-5:5))
		z = C.merge(a, b)
		n = length(z)
		# no truncation: exact
		f = C.Renorm(z, Val(n))
		@test exact(f) == exact(z)
		@test is_ulp_nonoverlapping(f)
		@test C.FastRenorm(z, Val(n)).val == f.val
		@test C.Renorm(z, C.CamparyFloat{n}(fill(T(7), n))).val == f.val
		# truncation to m < n terms: ulp-nonoverlapping, and the discarded tail of an
		# ulp-nonoverlapping expansion is at most ulp(f_{m-1}) / (1 - 2^(1-p))
		for m in 1:n-1
			g = C.Renorm(z, Val(m))
			@test g isa C.CamparyFloat{T, m}
			@test is_ulp_nonoverlapping(g)
			@test g.val == f.val[1:m]
			@test abs(exact(z) - exact(g)) <= BigFloat(eps(g[m])) / (1 - BigFloat(2)^(1 - prec(T)))
			@test C.FastRenorm(z, Val(m)).val == g.val
		end
	end
end

@testset "VecSumErr ($T)" for T in (Float64, Float32)
	for K in 1:maxlimbs(T), _ in 1:100
		x = rand_expansion(T, K)
		@test exact(C.VecSumErr(x, C.CamparyFloat{K}(zeros(T, K)))) == exact(x)
		y = C.CamparyFloat{K}(copy(x.val))
		@test exact(C.VecSumErr!(y)) == exact(x)
	end
end

@testset "merge ($T)" for T in (Float64, Float32)
	for _ in 1:200
		x, y = rand_expansion(T, rand(1:maxlimbs(T))), rand_expansion(T, rand(1:maxlimbs(T)); e = rand(-20:20))
		z = C.merge(x, y)
		@test length(z) == length(x) + length(y)
		@test issorted(abs.(z.val), rev = true)
		@test sort(z.val) == sort([x.val; y.val])
	end
end

# Theorem 3.3.2 of Popescu's thesis: |x + y - s| < (9/2) 2^(-(p-1) r) (|x| + |y|).
@testset "accurate addition, Algorithm 27 ($T)" for T in (Float64, Float32)
	p = prec(T)
	for r in 1:maxlimbs(T), _ in 1:300
		K, L = rand(1:maxlimbs(T)), rand(1:maxlimbs(T))
		x = rand_expansion(T, K)
		y = rand() < 0.2 ? C.CamparyFloat{K}(-x.val .* (1 + rand(T) * eps(T))) : rand_expansion(T, L; e = rand(-10:10))
		s = +(x, y, Val(r))
		@test s isa C.CamparyFloat{T, r}
		@test is_ulp_nonoverlapping(s)
		bound = BigFloat(9) / 2 * BigFloat(2)^(-(p - 1) * r) * (abs(exact(x)) + abs(exact(y)))
		@test abs(exact(x) + exact(y) - exact(s)) < bound
	end
	x, y = rand_expansion(T, 4), rand_expansion(T, 4)
	@test (x + y).val == (+(x, y, Val(4))).val
	@test (x + one(T)).val == (x + C.CamparyFloat{4}(one(T))).val
end

# Theorem 4.1.1 of Popescu's thesis, for Algorithm 40 (PAddition_safe).
@testset "data-parallel addition, Algorithm 40 ($T)" for T in (Float64, Float32)
	p = prec(T)
	for r in 1:maxlimbs(T), _ in 1:300
		x, y = rand_expansion(T, r), rand_expansion(T, r; e = rand(-3:3))
		s = C.PAddition(x, y)
		@test C._SIMDAddition(x, y).val == s.val
		@test s isa C.CamparyFloat{T, r}
		bound = BigFloat(2)^(-(p - 1) * r) * (abs(BigFloat(x[1])) + abs(BigFloat(y[1]))) *
			(BigFloat(2)^(2r) * (1 + BigFloat(2)^(-p)) + 1 / (1 - BigFloat(2)^(-(p - 1))))
		@test abs(exact(x) + exact(y) - exact(s)) < bound
	end
end

@testset "conversion ($T)" for T in (Float64, Float32)
	x = rand_expansion(T, 3)
	y = convert(C.CamparyFloat{T, 5}, x)
	@test y.val == [x.val; zeros(T, 2)]
	@test_throws ArgumentError convert(C.CamparyFloat{T, 2}, x)
	@test convert(C.CamparyFloat{T, 3}, T(2)).val == T[2, 0, 0]
end


@testset "is_ulp_nonoverlapping" begin
	@test is_ulp_nonoverlapping(C.CamparyFloat{3}([1.0, eps(1.0), eps(eps(1.0))]))
	@test !is_ulp_nonoverlapping(C.CamparyFloat{2}([1.0, nextfloat(eps(1.0))]))
	@test is_ulp_nonoverlapping(C.CamparyFloat{3}([1.0, 0.0, 0.0]))
	@test !is_ulp_nonoverlapping(C.CamparyFloat{3}([1.0, 0.0, 1e-300]))
end

# A literal scalar transcription of Algorithm 41 of Popescu's thesis (PMultiplication),
# lanes 0, ..., r-1 stored at indices 1, ..., r, line 18 read as x times y_{r-1}.
function alg41_scalar(x::Vector{T}, y::Vector{T}) where {T}
	r = length(x); s = zeros(T, r); π = zeros(T, r)
	for i in 0:r-2
		pe = [TwoProd(x[k], y[i+1]) for k in 1:r]; p = first.(pe); e = last.(pe)
		se = [TwoSum(s[k], p[k]) for k in 1:r]; s = first.(se); e′ = last.(se)
		π[i+1] = s[1]
		s = [s[2:end]; zero(T)]
		while any(!=(0), e)
			se = [TwoSum(s[k], e[k]) for k in 1:r]; s = first.(se); e = last.(se)
			e = [zero(T); e[1:end-1]]
		end
		while any(!=(0), e′)
			se = [TwoSum(s[k], e′[k]) for k in 1:r]; s = first.(se); e′ = last.(se)
			e′ = [zero(T); e′[1:end-1]]
		end
	end
	s = s .+ x .* y[r]
	π[r] = s[1]
	return π
end

@testset "data-parallel multiplication, Algorithm 41 ($T)" for T in (Float64, Float32)
	for r in 1:(T == Float64 ? 16 : 6), _ in 1:300
		x, y = rand_expansion(T, r), rand_expansion(T, r; e = rand(-5:5))
		π = C.PMultiplication(x, y)
		@test π isa C.CamparyFloat{T, r}
		@test π.val == alg41_scalar(x.val, y.val)
		r == 1 && @test π.val == [x[1] * y[1]]
	end
	# the cap on the inner loops: a non-finite input does not hang
	@test !all(isfinite, C.PMultiplication(C.CamparyFloat{3}(T[Inf, 0, 0]), C.CamparyFloat{3}(T[1, 0, 0])).val)
end

# Theorem 4.2.1 of Popescu's thesis does not hold for this r = 2 input; recorded so that
# a change to PMultiplication or to the reading of the theorem shows up here.
@testset "Theorem 4.2.1 counterexample" begin
	x = C.CamparyFloat{2}([1.2052776966458982, -1.980988149935599e-16])
	y = C.CamparyFloat{2}([0.016518450532383697, -3.2475787651990608e-18])
	@test is_ulp_nonoverlapping(x) && is_ulp_nonoverlapping(y)
	π = C.PMultiplication(x, y)
	p, r = big(53), big(2)
	a = abs(BigFloat(x[1]) * BigFloat(y[1]))
	rhs = a * BigFloat(2)^(-(p - 1) * r) * (BigFloat(128) / 127 * (2r - 1) - BigFloat(129) / 254 * r +
		BigFloat(2)^(-p - r + 2) * (r^2 - r) * factorial(r)^2)
	@test abs(exact(x) * exact(y) - exact(π)) > rhs
end

end
