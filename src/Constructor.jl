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



# Upward conversion of CamparyFloat types, i.e., to a bigger number of limbs: the extra
# limbs are zero, so the value is unchanged. Conversion to fewer limbs would round, and
# is refused; use `Renorm(x, Val(N))` for that.
function Base.convert(::Type{CamparyFloat{T, N}}, x::CamparyFloat{T, M}) where {T<:FloatTypes, M, N}
	N >= M || throw(ArgumentError("cannot convert a $M-limb CamparyFloat to $N limbs without rounding; use Renorm"))
	v = zeros(T, N)
	v[1:M] = x.val
	return CamparyFloat{N}(v)
end

function Base.convert(::Type{CuCamparyFloat{T, N}}, x::CamparyFloat{T, M}) where {T<:FloatTypes, M, N}
	N >= M || throw(ArgumentError("cannot convert a $M-limb CamparyFloat to $N limbs without rounding; use Renorm"))
	v = zeros(T, N)
	v[1:M] = x.val
	return CuCamparyFloat{N}(v)
end

"""
    is_ulp_nonoverlapping(x)

Whether the limbs of `x` satisfy `|x_{i+1}| <= ulp(x_i)` for every `i`, Definition 26 of
Boldo, Joldes, Muller, Popescu (ITP 2017); Definition 2.6 of Muller, Popescu, Tang
(ARITH 2016). The check is exact: `eps(z)` is `ulp(z)` for a normal `z`, and a zero limb
must be followed by zeros.
"""
function is_ulp_nonoverlapping(x::CamparyFloat{T,K}) where {T,K}
	for i in 1:K-1
		if x[i] == zero(T)
			x[i+1] == zero(T) || return false
		elseif !(abs(x[i+1]) <= eps(x[i]))
			return false
		end
	end
	return true
end
