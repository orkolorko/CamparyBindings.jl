# Data-parallel multiplication of floating-point expansions.
#
# [P17] V. Popescu, "Towards fast and certified multiple-precision libraries",
#       PhD thesis, ENS Lyon (2017), HAL tel-01534090, Section 4.2.
#
# Algorithm 41 of [P17] (PMultiplication), "r input - r output":
#
#  1: s <- (0, ..., 0)
#  2: pi <- (0, ..., 0)
#  3: for i <- 0 to r-2 do
#  4:    y' <- (y_i, ..., y_i)                    broadcast
#  5:    (p, e) <- 2Prod(x, y')
#  6:    (s, e') <- 2Sum(s, p)
#  7:    pi_i <- s_0
#  8:    s <- (s_1, ..., s_{r-1}, 0)               shift left
#  9:    while e != 0 do
# 10:       (s, e) <- 2Sum(s, e)
# 11:       e <- (0, e_0, ..., e_{r-2})            shift right
# 12:    end while
# 13:    while e' != 0 do
# 14:       (s, e') <- 2Sum(s, e')
# 15:       e' <- (0, e'_0, ..., e'_{r-2})         shift right
# 16:    end while
# 17: end for
# 18: p <- x . y
# 19: s <- s + p
# 20: pi_{r-1} <- s_0
#
# Line 18 is read as x times the broadcast y_{r-1}, which is what the CUDA code of the
# same algorithm in [P17], Figure 4.7, does (`y_i = shfl(y, R-1, R); p = x * y_i;`).
#
# The two while loops are run at most r times each. With finite s this changes nothing:
# lane 0 of e receives 0 at every shift, and 2Sum(s_0, 0) has zero error when s_0 is
# finite, so after t passes lanes 0, ..., t-1 of e are zero and the loop has ended by the
# r-th pass. The cap only matters when s holds an Inf or a NaN, where 2Sum(s_0, 0) returns
# a NaN error and the uncapped loop would not terminate.

using SIMD

# s <- (s_1, ..., s_{r-1}, 0) and e <- (0, e_0, ..., e_{r-2}), as SIMD shuffles.
@generated _shift_left(v::Vec{R,T}) where {R,T} =
	:(shufflevector(v, Vec{$R,$T}(zero($T)), Val($(Tuple(1:R)))))
@generated _shift_right(v::Vec{R,T}) where {R,T} =
	:(shufflevector(v, Vec{$R,$T}(zero($T)), Val($((R, (0:R-2)...)))))

@inline function _accumulate(s::Vec{R,T}, e::Vec{R,T}) where {R,T}
	z = Vec{R,T}(zero(T))
	for _ in 1:R
		any(e != z) || break
		s, e = TwoSum(s, e)
		e = _shift_right(e)
	end
	return s
end

"""
    PMultiplication(x, y)

Algorithm 41 of [P17] on two `R`-limb expansions, with the vector operations on
`SIMD.Vec{R}`: returns `R` limbs `pi` approximating `x * y`. The output limbs are not
renormalised and need not be ulp-nonoverlapping ([P17], p. 138). See the comment at the
top of `PMultiplication.jl` for the reading of line 18 and for the cap on the two inner
loops.

No error bound is attached. Theorem 4.2.1 of [P17] states, for ulp-nonoverlapping `x` and
`y`, `p >= 8` and no underflow or overflow,

    |xy - pi| <= |x_0 y_0| 2^(-(p-1) r) ((128/127)(2r - 1) - (129/254) r + 2^(-p-r+2) (r^2 - r) (r!)^2),

and the output of this routine, which agrees bit for bit with a literal scalar
transcription of Algorithm 41 (see the tests), exceeds it at `r = 2`: for
`x = (1.2052776966458982, -1.980988149935599e-16)` and
`y = (0.016518450532383697, -3.2475787651990608e-18)` in Float64,
`|xy - pi| = 2.2750 * 2^-104 |x_0 y_0|` against a right-hand side of
`2.0079 * 2^-104 |x_0 y_0|`. The proof charges the errors left out at the last step at
`Omega_r = 2^(-(r+1)p+2) |x_0 y_0| (r!)^2` each; in this example the rounding error of
`RN(x_1 y_0)` (left in lane 1), the 2Sum error of lane 0 (shifted to lane 1), and the
rounding error of line 19 are `0.14`, `0.39` and `0.78` times `2^-104 |x_0 y_0|`.
"""
function PMultiplication(x::CamparyFloat{T,R}, y::CamparyFloat{T,R}) where {T, R}
	xv = vload(Vec{R,T}, x.val, 1)
	s = Vec{R,T}(zero(T))
	π = zeros(T, R)
	for i in 1:R-1
		p, e = TwoProd(xv, Vec{R,T}(y[i]))
		s, e′ = TwoSum(s, p)
		π[i] = s[1]
		s = _shift_left(s)
		s = _accumulate(s, e)
		s = _accumulate(s, e′)
	end
	s = s + xv * Vec{R,T}(y[R])
	π[R] = s[1]
	return CamparyFloat{R}(π)
end
