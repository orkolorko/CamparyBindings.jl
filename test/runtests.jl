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

# Proposition 34 of Boldo-Joldes-Muller-Popescu (ITP 2017): |f_{j+1}| <= ulp(f_j).
is_ulp_nonoverlapping(f) = all(abs(f[j+1]) <= eps(f[j]) for j in 1:length(f)-1)

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

end
