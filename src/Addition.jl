# Addition of floating-point expansions, "accurate" variant.
#
# [P17] V. Popescu, "Towards fast and certified multiple-precision libraries",
#       PhD thesis, ENS Lyon (2017), HAL tel-01534090.

"""
    merge(x, y)
    merge(x, y, Val(n0), Val(m0))

Algorithm 28 of [P17]: merges the limbs of `x` and `y` into one array in decreasing order
of magnitude, taking `y_j` when `|y_j| > |x_i|` and `x_i` otherwise. The second form
merges only the first `n0` limbs of `x` and the first `m0` limbs of `y`.
"""
merge(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}) where {T,K,L} = merge(x, y, Val(K), Val(L))

function merge(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}, ::Val{N}, ::Val{M}) where {T,K,L,N,M}
	@assert N <= K && M <= L
	f = Vector{T}(undef, N+M)
	i, j = 1, 1
	for t in 1:N+M
		if i > N || (j <= M && abs(y[j]) > abs(x[i]))
			f[t] = y[j]; j += 1
		else
			f[t] = x[i]; i += 1
		end
	end
	return CamparyFloat{N+M}(f)
end

"""
    +(x, y, Val(r))

Algorithm 27 of [P17] (`Addition_accurate`): truncates `x` and `y` to at most `r` limbs,
merges them with [`merge`](@ref) and renormalises to `r` limbs with `Renorm`. For
ulp-nonoverlapping `x` and `y` and `p >= 4`, without underflow or overflow, the result
`s` is ulp-nonoverlapping and `|x + y - s| < (9/2) 2^(-(p-1) r) (|x| + |y|)`
(Theorem 3.3.2 of [P17]).

`x + y` with both of `K` limbs is `+(x, y, Val(K))`.
"""
Base.:+(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}, ::Val{R}) where {T,K,L,R} =
	Renorm(merge(x, y, Val(min(K, R)), Val(min(L, R))), Val(R))
Base.:+(x::CamparyFloat{T,K}, y::CamparyFloat{T,K}) where {T,K} = Base.:+(x, y, Val(K))
