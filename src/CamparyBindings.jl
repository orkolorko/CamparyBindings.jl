module CamparyBindings
using StaticArrays

export CamparyFloat, FastTwoSum, FastTwoDiff, TwoSum, TwoProd, VecSum, FastVecSum, Sum

struct CamparyFloat{T, N} <: Real
	val::SVector{N, T}
	function CamparyFloat{N}(v::V) where {T, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T, N}(SVector{N, T}(v))
	end
end

Base.getindex(x::CamparyFloat{T,K}, i::Integer) where {T, K} = getindex((x.val), i)
Base.length(x::CamparyFloat{T,K}) where {T, K} = length(x.val)

@inline FastTwoSum(x, y) = (x+y, y-((x+y)-x))
@inline FastTwoDiff(x,y) = (x-y, (x-(x-y))-y)

@inline function TwoSum(x, y) 
  	s = x+y
  	aa = s-y
  	bb = s-aa
  	da = x-aa
  	db = y-bb
  	return s, da+db
end

@inline function TwoDiff(a, b)
	s = a-b
	bb = s-a
	aa = s-bb
	da = a-aa
	db = b+bb
	return s, da-db
end

@inline TwoProd(a, b) = (a*b, fma(a, b, -a*b))

merge(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}) where {T,K,L} = CamparyFloat{K+L}([x.val ; y.val]) 


function VecSum(x::CamparyFloat{T,K}) where {T, K}
	v = Array{T}(undef, K)
	v[K] = x[K]
	for i in 1:K-1
		@inbounds v[K-i], v[K-i+1] = TwoSum(x[K-i], v[K-i+1])
	end
	return CamparyFloat{K}(v)
end

function FastVecSum(x::CamparyFloat{T,K}) where {T, K}
	v = Array{T}(undef, K)
	v[K] = x[K]
	for i in 1:K-1
		@inbounds v[K-i], v[K-i+1] = FastTwoSum(x[K-i], v[K-i+1])
	end
	return CamparyFloat{K}(v)
end

struct Sum{L}
end

function (p::Sum{L})(x) where {L}
	return x
end

#function (p::Sum{L})(x::CamparyFloat{T, K}) where {T,K,L}
#	res = Array{T}(undef, L)
	#for i in 1:K-L
	#	x = VecSum(x)
	#end
	#for k in 0:L-2
	#	x[1:K-k] = VecSum(x[1:K-k])
	#	res[k+1] = x[K-k]
	#end
	#res[L] = sum(x[1:K-L+1])
#	return res
#end



### Is it worth to use generated????
#@generated function VecSum(x::CamparyFloat{T,K}) where {T, K}
#	ex = quote end
#	push!(ex.args, quote
#			v = Array{$T}(undef, $K)
#			v[$K] = x[$K]
#		end)
#	for i in K-1:-1:1
#		push!(ex.args, 
#			quote
#				v[$i], v[$(i+1)] = TwoSum(x[$i], v[$(i+1)])
#			end)
#	end
#	push!(ex.args, 
#		quote
#			return CamparyFloat{$K}(v)
#		end)
#	return ex
#end

#@generated function FastVecSum(x::CamparyFloat{T,K}) where {T, K}
#	ex = quote end
#	push!(ex.args,
#		quote
#			v = Array{$T}(undef, $K)
#			v[$K] = x[$K]
#		end
#		)
#
#	for i in K-1:-1:1
#		push!(ex.args, 
#			quote
#				v[$i], v[$(i+1)] = FastTwoSum(x[$i], v[$(i+1)])
#			end)
#	end
#	push!(ex.args, 
#		quote
#			return CamparyFloat{$K}(v)
#		end)
#	return ex
#end


end # module
