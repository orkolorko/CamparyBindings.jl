FloatTypes = Union{Float32, Float64}

struct CamparyFloat{T, N} <: Real
	val::Vector{T}
	function CamparyFloat{N}(v::V) where {T<:FloatTypes, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T,N}(v)
	end
end

#Static Version
struct SCamparyFloat{T, N} <: Real
	val::SVector{N, T}
	function SCamparyFloat{N}(v::V) where {T<:FloatTypes, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T, N}(SVector{N, T}(v))
	end
end

CamparyFloat{N}(x::T) where {T<:FloatTypes,N} = (v=zeros(T, N); v[1]=x; return CamparyFloat{N}(v))

CamparyFloatTypes{T, N} = Union{CamparyFloat{T, N}, SCamparyFloat{T, N}}

Base.getindex(x::CamparyFloatTypes{T, K}, i) where {T, K} = CamparyFloat{length(i)}(getindex((x.val), i))
Base.getindex(x::CamparyFloatTypes{T, K}, i::Integer) where {T, K} = getindex((x.val), i)


Base.length(x::CamparyFloatTypes{T, K}) where {T, K} = length(x.val)

### Maybe introduce some sanity checks???
Base.setindex!(x::CamparyFloat{T, K}, y::T, i) where {T, K} = setindex!(x.val, y, i) 
Base.setindex!(x::CamparyFloat{T, K}, v::Vector, i) where {T, K} = setindex!(x.val, v, i) 
Base.setindex!(x::CamparyFloat{T, K}, y::CamparyFloatTypes{T, M}, i) where {T, K, M} = setindex!(x.val, y.val, i) 

Base.convert(::Type{CamparyFloat{T₁, N}}, x::T₂) where {T₁<:FloatTypes, N, T₂<:Real} = CamparyFloat{N}(convert(T₁, x))
# we allow upward conversion of CamparyFloat types, i.e., bigger number of limbs
function Base.convert(::Type{CamparyFloat{T, N}}, x::T₂) where {T<:Union{Float32, Float64}, M, N, T₂<:CamparyFloat{T,M}}
	if N>=M
		if M == 1
			exp = quote
				v = zeros($T, $N)
				v[1] = x[1]
				return CamparyFloat{$N}(v)
			end
			return exp
		end
		if M != 1
			exp = quote
				v = zeros($T, $N)
				v[1:$M] = x[1:$M]
				return CamparyFloat{$N}(v)
			end
			return exp
		end
	else
		@error "Invalid conversion since $M > $N"
	end
end



