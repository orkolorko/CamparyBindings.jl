# Data-parallel addition of floating-point expansions.
#
# [P17] V. Popescu, "Towards fast and certified multiple-precision libraries",
#       PhD thesis, ENS Lyon (2017), HAL tel-01534090.
#
# Both routines below are Algorithm 40 of [P17] (`PAddition_safe`), with the first loop
# read as i = 1, ..., r-1: Algorithm 40 prints "for i <- 1 to r", but x_r and y_r do not
# exist, and the text says that after the inputs are consumed the errors are propagated
# "for another r - 1 steps", which is the r - 2 steps of the second loop plus the final
# plain addition. Theorem 4.1.1 of [P17]: for ulp-nonoverlapping x and y, without
# underflow or overflow,
#     |x + y - s| < 2^(-(p-1) r) (|x_0| + |y_0|) (2^(2r) (1 + 2^-p) + 1/(1 - 2^(-(p-1)))).

using SIMD

function VecTwoSum!(v, w, ::Val{R}, s, e) where {R}
	@inbounds for i in 1:R
		s[i], e[i] = TwoSum(v[i], w[i])
	end
end

"""
    PAddition(x, y)

Algorithm 40 of [P17] on `R`-limb expansions, with the vector operations written as
loops over the limbs. See the comment at the top of `PAddition.jl` for the loop bound and
the error bound.
"""
function PAddition(x::CamparyFloat{T,R}, y::CamparyFloat{T,R}) where {T, R}
	s = zeros(T, R)
	s[1] = x[1]
	e = zeros(T, R)
	e[1] = y[1]

	VecTwoSum!(s, e, Val(R), s, e)

	for i in 2:R
		e′ = [x[i] ; e[1:end-1]]
		VecTwoSum!(s, e′, Val(R), s, e)
		e′ = [y[i] ; e[1:end-1]]
		VecTwoSum!(s, e′, Val(R), s, e)
	end

	for i in 1:R-2
		e′ = [zero(T) ; e[1:end-1]]
		VecTwoSum!(s, e′, Val(R), s, e)
	end
	e′ = [zero(T) ; e[1:end-1]]
	s += e′
	return CamparyFloat{R}(s)
end

"""
    _SIMDAddition(x, y)

Algorithm 40 of [P17] on `R`-limb expansions with SIMD.jl vectors: the shift right and
insert of lines 5, 7, 11 and 14 is a `shufflevector`. Gives the same result as
[`PAddition`](@ref) bit for bit.
"""
@generated function _SIMDAddition(x::CamparyFloat{T,R}, y::CamparyFloat{T,R}) where {T, R}
	expr=
	quote
	s = Vec{$R, $T}(zero($T))
	e = Vec{$R, $T}(zero($T))
	z = Vec{$R, $T}(zero($T)) # the zero vector
	x_ass = vload(Vec{$R,$T}, x.val, 1) #loads the values into the processor
	y_ass = vload(Vec{$R,$T}, y.val, 1)
	end

	shift = NTuple{R-1,Int}(R+i for i in 0:R-2)
	mask = (0, shift...) #creates the tuple

	push!(expr.args,
	quote
	s = shufflevector(x_ass, s, Val($mask))
	e = shufflevector(y_ass, e, Val($mask))
	s, e = TwoSum(s, e)
	end)

	for i in 1:R-1
		mask = (i, shift...) #creates the tuple
		push!(expr.args,
		quote
		e′ = shufflevector(x_ass, e, Val($mask))
		s, e = TwoSum(s, e′)
		e′ = shufflevector(y_ass, e, Val($mask))
		s, e = TwoSum(s, e′)
		end)
	end

	mask = (0, shift...) #creates the tuple
	for i in 1:R-2
		push!(expr.args,
		quote
		e′ = shufflevector(z, e, Val($mask))
		s,e = TwoSum(s, e′)
		end
		)
	end
	push!(expr.args,
	quote
	e′ = shufflevector(z, e, Val($mask))
	s += e′
	res = zeros($T, $R)
	vstore(s, res, 1)
	return CamparyFloat{$R}(res)
	end)
	return expr
end
