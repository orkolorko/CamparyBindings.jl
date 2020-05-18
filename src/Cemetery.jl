#Static Version
struct SCamparyFloat{T, N} <: Real
	val::SVector{N, T}
	function SCamparyFloat{N}(v::V) where {T<:FloatTypes, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T, N}(SVector{N, T}(v))
	end
end



# This file contains code that may be worth preserving as a version check

function old_VecSum!(x::CamparyFloatTypes{T, K}) where {T, K}
	s = x[K]
	for i in 1:K-1
		s, x[K-i+1] = TwoSum(x[K-i], s)
	end
	x[1] = s
	return x
end

function old_FastVecSum!(x::CamparyFloatTypes{T, K}) where {T, K}
	s = x[K]
	for i in 1:K-1
		s, x[K-i+1] = FastTwoSum(x[K-i], s)
	end
	x[1] = s
	return x
end

function old_VecSum(x::CamparyFloatTypes{T, K}) where {T, K}
	v = Array{T}(undef, K)
	v[K] = x[K]
	for i in 1:K-1
		@inbounds v[K-i], v[K-i+1] = TwoSum(x[K-i], v[K-i+1])
	end
	return CamparyFloat{K}(v)
end

function old_FastVecSum(x::CamparyFloatTypes{T, K}) where {T, K}
	v = Array{T}(undef, K)
	v[K] = x[K]
	for i in 1:K-1
		@inbounds v[K-i], v[K-i+1] = FastTwoSum(x[K-i], v[K-i+1])
	end
	return CamparyFloat{K}(v)
end

function old_Sum(x::CamparyFloat{T, K}, ::Val{L}) where {T,K,L}
	res = Array{T}(undef, L)
	for i in 1:K-L
		VecSum!(x)
	end
	for k in 0:L-2
		x[k+1:K] = VecSum!(x[k+1:K])
		res[k+1] = x[k]
	end
	res[1] = reduce(+, (x.val)[1:K-L+1])
	return res
end

old_Sum1(x::CamparyFloat{T, K}) where {T,K} = Sum(x, Val(1))
old_SumK(x::CamparyFloat{T, K}) where {T,K} = Sum(x, Val(K))

islessequal(::Val{K}, ::Val{L}) where {K, L} = Val(K<=L)

function FastVecSumErrBranch!(x::CamparyFloat{T,K}, ::Val{R}) where {T, K, R}
 	return _FastVecSumErrBranch!(x, Val(R), islessequal(Val(R), Val(K)))
end

function _FastVecSumErrBranch!(x::CamparyFloat{T,K}, ::Val{R}, ::Val{true}) where {T, K, R}
 	ptr = 1
 	i = 2
 	e = x[1]

 	while (i<=K && ptr<=R) 
 	 	x[ptr], e = FastTwoSum(e, x[i])
 	 	i=i+1
		
 	 	if e == zero(T)
 	 		e = x[ptr]
 	 	else
 	 		ptr+= 1
 	 	end
 	end
 	if ptr<=R &&  e!=zero(T)
 	 	x[ptr] = e
 	 	ptr+=1
 	end
 	for i in ptr:R
 	 	x[i] = zero(T)
 	end
 	return x
end

function _FastVecSumErrBranch!(x::CamparyFloat{T,K}, ::Val{R}, ::Val{false}) where {T, K, R}
	@error "The result size is too big, $K < $R"
end

## Is it worth to use generated????
@generated function VecSum(x::CamparyFloat{T,K}) where {T, K}
	ex = quote end
	push!(ex.args, quote
			v = Array{$T}(undef, $K)
			v[$K] = x[$K]
		end)
	for i in K-1:-1:1
		push!(ex.args, 
			quote
				v[$i], v[$(i+1)] = TwoSum(x[$i], v[$(i+1)])
			end)
	end
	push!(ex.args, 
		quote
			return CamparyFloat{$K}(v)
		end)
	return ex
end

@generated function FastVecSum(x::CamparyFloat{T,K}) where {T, K}
	ex = quote end
	push!(ex.args,
		quote
			v = Array{$T}(undef, $K)
			v[$K] = x[$K]
		end
		)

	for i in K-1:-1:1
		push!(ex.args, 
			quote
				v[$i], v[$(i+1)] = FastTwoSum(x[$i], v[$(i+1)])
			end)
	end
	push!(ex.args, 
		quote
			return CamparyFloat{$K}(v)
		end)
	return ex
end
