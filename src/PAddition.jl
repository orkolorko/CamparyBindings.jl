using SIMD


unzip(v)= ([val[1] for val in v], [val[2] for val in v])

function VecTwoSum!(v, w, ::Val{R}, s, e) where {R}
	@inbounds for i in 1:R
		s[i], e[i] = TwoSum(v[i], w[i])
	end
end


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
		e′ = [0.0 ; e[1:end-1]]
		VecTwoSum!(s, e′, Val(R), s, e)
	end
	e′ = [0.0 ; e[1:end-1]]
	s += e′
	return CamparyFloat{R}(s)
end

function SIMDTwoSum!(v, w, ::Val{R}, s, e) where {T, R}
	s = v+w
	aa = s-
  		bb = s[lane+i]-aa
  		da = v[lane+i]-aa
  		db = w[lane+i]-bb
  		e[lane+i] = da+db
end


is_even(::Val{R}) where {R} = Val(R%2==0)

@generated function _SIMDAddition(x::CamparyFloat{T,R}, y::CamparyFloat{T,R}, ::Val{true}) where {T, R}
	expr=
	quote
	s = Vec{$R, $T}(0.0)
	e = Vec{$R, $T}(0.0)
	z = Vec{$R, $T}(0.0) # the zero vector
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
	res = zeros(Float64, R)
	vstore(s, res, 1)
	return CamparyFloat{$R}(res)
	end)
	return expr
end

