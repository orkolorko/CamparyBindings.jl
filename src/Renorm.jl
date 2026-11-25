# Renormalisation of floating-point expansions.
#
# [BJMP17] S. Boldo, M. Joldes, J.-M. Muller, V. Popescu, "Formal verification of a
#          floating-point expansion renormalization algorithm", Interactive Theorem
#          Proving (ITP 2017), LNCS, pp. 98--113, doi:10.1007/978-3-319-66107-0_7.
# [JMMP16] M. Joldes, O. Marty, J.-M. Muller, V. Popescu, "Arithmetic algorithms for
#          extended precision using floating-point expansions", IEEE Transactions on
#          Computers 65 (2016), pp. 1197--1210, doi:10.1109/tc.2015.2441714.
# [P17]    V. Popescu, "Towards fast and certified multiple-precision libraries",
#          PhD thesis, ENS Lyon (2017), HAL tel-01534090.

"""
    FastVecSum(x)
    FastVecSum(x, r)
    FastVecSum!(x)

[`VecSum`](@ref) with `FastTwoSum` in place of `TwoSum`. [BJMP17] proves the substitution
admissible for input that overlaps by at most `d <= p-2` digits.
"""
function FastVecSum(x::CamparyFloatTypes{T, K}, r::CamparyFloat{T, K}) where {T, K}
	s = x[K]
	for i in K-1:-1:1
		s, r[i+1] = FastTwoSum(x[i], s)
	end
	r[1] = s
	return r
end
FastVecSum(x::CamparyFloatTypes{T, K}) where {T, K} = FastVecSum(x, CamparyFloat{K}(Vector{T}(undef, K)))
FastVecSum!(x::CamparyFloat{T, K}) where {T, K} = FastVecSum(x, x)

"""
    VecSum(x)
    VecSum(x, r)
    VecSum!(x)

Algorithm 4 of [BJMP17] (Algorithm 21 of [P17]): the chain of `TwoSum` from the least
significant term, `s_{n-1} = x_{n-1}`, `(s_i, e_{i+1}) = 2Sum(x_i, s_{i+1})`, `e_0 = s_0`.
Error-free: the exact sum of the output equals the exact sum of `x`. The second form
writes into `r`, the third into `x`.
"""
function VecSum(x::CamparyFloatTypes{T, K}, r::CamparyFloat{T, K}) where {T, K}
	s = x[K]
	for i in K-1:-1:1
		s, r[i+1] = TwoSum(x[i], s)
	end
	r[1] = s
	return r
end
VecSum(x::CamparyFloatTypes{T, K}) where {T, K} = VecSum(x, CamparyFloat{K}(Vector{T}(undef, K)))
VecSum!(x::CamparyFloat{T, K}) where {T, K} = VecSum(x, x)


# This functions add a CamparyFloat and a float and return a CamparyFloat{K+1}

# We use the following functions to forbid the compilation of the function
# if r is not one limb bigger than x
issuccessor(::Val{K}, ::Val{L}) where {K, L} = Val(K-1==L)

function VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}) where {T,K,L}
	return _VecSum4Add1(x, y, r, issuccessor(Val(L), Val(K)))
end

function _VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}, ::Val{true}) where {T,K, L}
	s = y
	for i in K:-1:1
		s, r[i+1] = TwoSum(x[i], s)
	end
	r[1] = s
	return r
end

function _VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}, ::Val{false}) where {T,K, L}
	throw(DimensionMismatch("the output has $L limbs, it must have $K + 1"))
end

VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T) where {T, K} =
	_VecSum4Add1(x, y, CamparyFloat{K+1}(Vector{T}(undef, K+1)), Val(true))

"""
    VecSumErrBranch(e, Val(m))
    VecSumErrBranch(e, r)

Algorithm 5 of [BJMP17] (Algorithm 22 of [P17]): a chain of `TwoSum` from the most
significant term that propagates the error, or the sum when the error is zero, and stops
once `m` output terms are filled. For S-nonoverlapping input the output `f` satisfies
`|f_{j+1}| <= ulp(f_j)` (Proposition 34 of [BJMP17]). Unfilled output limbs are zero. The
second form writes into `r` and takes `m` from its number of limbs.
"""
VecSumErrBranch(x::CamparyFloat{T,K}, ::Val{R}) where {T, K, R} = VecSumErrBranch(x, CamparyFloat{R}(zeros(T, R)))

function VecSumErrBranch(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T, K, R}
	fill!(r.val, zero(T))
	j = 1
	err = x[1]
	for i in 1:K-1
		r[j], err = TwoSum(err, x[i+1])
		if err != zero(T)
			if j>=R
				return r
			end
			j+= 1
		else
			err = r[j]
		end
	end
	if err !=zero(T) && j<=R
		r[j]=err
	end
	return r
end

"""
    FastVecSumErrBranch(e, Val(m))
    FastVecSumErrBranch(e, r)

`VecSumErrBranch` with `FastTwoSum` in place of `TwoSum`, which [BJMP17] proves
admissible for this level. The control flow is that of Algorithm 5 written as a single
loop; unfilled output limbs are zero. (Author's note: even if it is called Fast, the speed
seems to be the same.)
"""
function FastVecSumErrBranch(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T, K, R}
 	# in this case the check is done in the calling function

 	ptr = 1
 	i = 2
 	e = x[1]

 	while i<=K && ptr<=R
 		r[ptr], e = FastTwoSum(e, x[i])
 		i+=1
		
 		if e == zero(T)
 			e = r[ptr]
 		else
 			ptr+=1
 		end
 	end
 	if ptr<=R &&  e!=zero(T)
 		r[ptr] = e
 		ptr+=1
 	end
 	for i in ptr:R
 		r[i] = zero(T)
 	end
 	return r
end

function FastVecSumErrBranch(x::CamparyFloat{T,K}, ::Val{R}) where {T, K, R}
 	# no check needed, the allocation is here
 	r = Vector{T}(undef, R)	

 	ptr = 1
 	i = 2
 	e = x[1]

 	while i<=K && ptr<=R
 		r[ptr], e = FastTwoSum(e, x[i])
 		i+= 1
		
 		if e == zero(T)
 			e = r[ptr]
 		else
 			ptr+=1
 		end
 	end
 	if ptr<=R &&  e!=zero(T)
 		r[ptr] = e
 		ptr+=1
 	end
	for i in ptr:R
 		r[i] = zero(T)
 	end
 	return CamparyFloat{R}(r)
end

"""
    VecSumErr!(f)
    VecSumErr(f, g)

Algorithm 8 of [JMMP16], the third level of the renormalisation: `r_0 = f_0`,
`(g_i, r_{i+1}) = 2Sum(r_i, f_{i+1})` for `i = 0, ..., m-1`, `g_m = r_m`. Error-free.
"""
function VecSumErr(x::CamparyFloat{T,K}, r::CamparyFloat{T,K}) where {T,K}
	err = x[1]
	for i in 1:K-1
		r[i], err = TwoSum(err, x[i+1])
	end
	r[K] = err
	return r
end
VecSumErr!(x::CamparyFloat{T,K}) where {T,K} = VecSumErr(x, x)

"""
    Renorm(x, Val(m))
    Renorm(x, r)

Algorithm 3 of [BJMP17] (Algorithm 24 of [P17]): `VecSum` followed by `VecSumErrBranch`.
If the terms of `x` overlap by at most `d <= p-2` digits in Priest's sense (which
requires them in decreasing order of magnitude, up to interleaving zeros), the result is
a truncation to `m` terms of an ulp-nonoverlapping expansion with the same exact sum as
`x` (Proposition 32 of [BJMP17], with a Coq proof there).
"""
Renorm(x::CamparyFloat{T,K}, ::Val{R}) where {T,K,R} = VecSumErrBranch(VecSum(x), Val(R))
Renorm(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T,K,R} = VecSumErrBranch(VecSum(x), r)

"""
    FastRenorm(x, Val(m))
    FastRenorm(x, r)

`Renorm` with `FastTwoSum` throughout, which [BJMP17] proves admissible under the same
hypothesis.
"""
FastRenorm(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T,K,R} = FastVecSumErrBranch(FastVecSum(x), r)
FastRenorm(x::CamparyFloat{T,K}, ::Val{R}) where {T,K,R} = FastVecSumErrBranch(FastVecSum(x), Val(R))