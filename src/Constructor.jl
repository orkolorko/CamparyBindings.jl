FloatTypes = Union{Float32, Float64}

struct CamparyFloat{T, N} <: Real
	val::Vector{T}
	function CamparyFloat{N}(v::V) where {T<:FloatTypes, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T,N}(v)
	end
end

struct CuCamparyFloat{T, N} <: Real
	val::CuVector{T}
	function CuCamparyFloat{N}(v::V) where {T<:FloatTypes, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T,N}(v)
	end
end



CamparyFloat{N}(x::T) where {T<:FloatTypes,N} = (v=zeros(T, N); v[1]=x; return CamparyFloat{N}(v))
CuCamparyFloat{N}(x::T) where {T<:FloatTypes,N} = (v=zeros(T, N); v[1]=x; return CuCamparyFloat{N}(v))

CamparyFloatTypes{T, N} = Union{CamparyFloat{T, N}, CuCamparyFloat{T, N}}

Base.getindex(x::CamparyFloatTypes{T, K}, i) where {T, K} = getindex((x.val), i)
Base.getindex(x::CamparyFloatTypes{T, K}, i::Integer) where {T, K} = getindex((x.val), i)

Base.length(x::CamparyFloatTypes{T, K}) where {T, K} = length(x.val)

### Maybe introduce some sanity checks???
Base.setindex!(x::CamparyFloat{T, K}, y::T, i) where {T, K} = setindex!(x.val, y, i) 
Base.setindex!(x::CamparyFloat{T, K}, v::Vector, i) where {T, K} = setindex!(x.val, v, i) 
Base.setindex!(x::CamparyFloatTypes{T, K}, y::CamparyFloatTypes{T, M}, i) where {T, K, M} = setindex!(x.val, y.val, i) 

Base.convert(::Type{CamparyFloat{T₁, N}}, x::T₂) where {T₁<:FloatTypes, N, T₂<:Real} = CamparyFloat{N}(convert(T₁, x))
Base.convert(::Type{CuCamparyFloat{T₁, N}}, x::T₂) where {T₁<:FloatTypes, N, T₂<:Real} = CuCamparyFloat{N}(convert(T₁, x))
Base.promote_rule(::Type{CamparyFloat{T₁, N}}, ::Type{T₂}) where {T₁<:FloatTypes, N, T₂<:Real} =
CamparyFloat{T₁,N}
Base.promote_rule(::Type{CuCamparyFloat{T₁, N}}, ::Type{T₂}) where {T₁<:FloatTypes, N, T₂<:Real} =
CuCamparyFloat{T₁,N}



# we allow upward conversion of CamparyFloat types, i.e., bigger number of limbs
@generated function Base.convert(::Type{CamparyFloat{T, N}}, x::T₂) where {T<:Union{Float32, Float64}, M, N, T₂<:CamparyFloat{T,M}}
	if N>=M
		if M == 1
			exp = quote
				v = zeros($T, $N)
				v[1] = x.val[1]
				return CamparyFloat{$N}(v)
			end
			return exp
		end
		if M != 1
			exp = quote
				v = zeros($T, $N)
				v[1:$M] = x.val[1:$M]
				return CamparyFloat{$N}(v)
			end
			return exp
		end
	else
		@error "Invalid conversion since $M > $N"
	end
end

@generated function Base.convert(::Type{CuCamparyFloat{T, N}}, x::T₂) where {T<:Union{Float32, Float64}, M, N, T₂<:CamparyFloat{T,M}}
	if N>=M
		if M == 1
			exp = quote
				v = zeros($T, $N)
				v[1] = x.val[1]
				return CuCamparyFloat{$N}(v)
			end
			return exp
		end
		if M != 1
			exp = quote
				v = zeros($T, $N)
				v[1:$M] = x.val[1:$M]
				return CuCamparyFloat{$N}(v)
			end
			return exp
		end
	else
		@error "Invalid conversion since $M > $N"
	end
end


