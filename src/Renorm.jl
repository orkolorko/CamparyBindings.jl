function FastVecSum!(x::CamparyFloatTypes{T, K}) where {T, K}
	s, x[K] = FastTwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, x[i+1] =FastTwoSum(x[i], s)
	end
	x[1], x[2] = FastTwoSum(x[1], s)
	return x
end

function FastVecSum(x::CamparyFloatTypes{T, K}) where {T, K}
	r = Vector{T}(undef, K)
	s, r[K] = FastTwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, r[i+1] =FastTwoSum(x[i], s)
	end
	r[1], r[2] = FastTwoSum(x[1], s)
	return CamparyFloat{K}(r)
end

function FastVecSum(x::CamparyFloatTypes{T, K}, r::CamparyFloat{T, K}) where {T, K}
	s, r[K] = FastTwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, r[i+1] =FastTwoSum(x[i], s)
	end
	r[1], r[2] = FastTwoSum(x[1], s)
	return r
end

function VecSum!(x::CamparyFloatTypes{T, K}) where {T, K}
	s, x[K] = TwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, x[i+1] =TwoSum(x[i], s)
	end
	x[1], x[2] = TwoSum(x[1], s)
	return x
end

function VecSum(x::CamparyFloatTypes{T, K}) where {T, K}
	r = Vector{T}(undef, K)
	s, r[K] = TwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, r[i+1] =TwoSum(x[i], s)
	end
	r[1], r[2] = TwoSum(x[1], s)
	return CamparyFloat{K}(r) 
end

function VecSum(x::CamparyFloatTypes{T, K}, r::CamparyFloat{T, K}) where {T, K}
	s, r[K] = TwoSum(x[K-1], x[K])
	for i in K-2:-1:2
		s, r[i+1] =TwoSum(x[i], s)
	end
	r[1], r[2] = TwoSum(x[1], s)
	return r
end


# This functions add a CamparyFloat and a float and return a CamparyFloat{K+1}

# We use the following functions to forbid the compilation of the function
# if r is not one limb bigger than x
issuccessor(::Val{K}, ::Val{L}) where {K, L} = Val(K-1==L)

function VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}) where {T,K,L}
	return _VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}, issuccessor(Val(L), Val(K)))
end

function _VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}, ::Val{true}) where {T,K, L}
	s, r[K+1] = TwoSum(x[K], y)
	for i in K-1:-1:2
		s, r[i+1] = TwoSum(x[i], s)
	end
	r[1], r[2] = TwoSum(x[1], s)
	return r
end 

function _VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T, r::CamparyFloatTypes{T, L}, ::Val{false}) where {T,K, L}
	@error "Incompatible sizes $L != $K +1" 
end 


function VecSum4Add1(x::CamparyFloatTypes{T, K}, y::T) where{T, K}
	r = Vector{T}(undef, K+1)
	s, r[K+1] = TwoSum(x[K], y)
	for i in K-1:-1:2
		s, r[i+1] = TwoSum(x[i], s)
	end
	r[1], r[2] = TwoSum(x[1], s)
	return CamparyFloat{K+1}(r)
end 

# this is another form of the same algorithm, from 
# Formal Verification of a Floating-Point Expansion
# Renormalization Algorithm by 
# Sylvie Boldo, Mioara Joldes, Jean-Michel Muller, Valentina Popescu
#
## I could reuse some code here
function VecSumErrBranch(x::CamparyFloat{T,K}, ::Val{R}) where {T, K, R}
	# no check needed, the allocation is here
 	r = zeros(T, R)	
	
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
	return CamparyFloat{R}(r)
end

function VecSumErrBranch(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T, K, R}
	# no check needed, the allocation is here
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

### Even if it is called Fast, the speed seems to be the same
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

function VecSumErr!(x::CamparyFloat{T,K}) where{T,K}
	for i in 2:K-1
		x[i-1], err = TwoSum(err, x[i])
	end
	x[K-1], x[K] = TwoSum(err, x[K])
	return x 
end

function VecSumErr(x::CamparyFloat{T,K}, r::CamparyFloat{T,K}) where{T,K}
	for i in 2:K-1
		r[i-1], err = TwoSum(err, x[i])
	end
	r[K-1], r[K] = TwoSum(err, x[K])
	return r
end

# Renormalization algorithm, from 
# Formal Verification of a Floating-Point Expansion
# Renormalization Algorithm by 
# Sylvie Boldo, Mioara Joldes, Jean-Michel Muller, Valentina Popescu

function Renorm(x::CamparyFloat{T,K}, ::Val{R}) where {T,K,R}
	x = VecSum(x)
	return  VecSumErrBranch(x, Val(R))
end

function Renorm(x::CamparyFloat{T,K}, r::CamparyFloat{T,R}) where {T,K,R}
	x = VecSum(x)
	return  VecSumErrBranch(x, r)
end
