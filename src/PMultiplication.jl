function _SIMDMultiplication(x::CamparyFloat{T,R}, y::CamparyFloat{T,R}, ::Val{true}) where {T, R}
	s = Vec{R, T}(0.0)
	z = Vec{R, T}(0.0) # the zero vector
	x_ass = vload(Vec{$R,$T}, x.val, 1) #loads the values into the processor
	y_ass = vload(Vec{$R,$T}, y.val, 1)

	π = zeros(T, R) #vector storing the result

	shift_sl = NTuple{R-1,Int}(i for i in 1:R-1)
	mask_sl = (shift_sl..., R)
	
	shift_sr = NTuple{R-1,Int}(i for i in 0:R-2)	
	mask_sr = (R, shift_sr...)

	for i in 0:R-1
		y′ = Vec{R,T}(y[i])
		p, e = TwoProd(x_ass, y')
		s, e′ = TwoSum(s, p)
		π[i+1] = s[0]
		s = shufflevector(s, e, Val(mask))
	end

	while !all(e == z)
		s, e = TwoSum(s, e)
		e = shufflevector(e, z, mask_sr)
	end
	while !all(e′ == z)
		s, e′ = TwoSum(s, e′)
		e′ = shufflevector(e′, z, mask_sr)
	end
	p = x_ass*y_ass
	s = s+p
	π[R] = s[0]
	return CamparyFloat{R}(π)
end
