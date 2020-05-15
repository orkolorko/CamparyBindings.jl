module CamparyBindings
using StaticArrays

export CamparyFloat, merge

struct CamparyFloat{T,N}
	val::SVector{N, T}
	function CamparyFloat{N}(v::V) where {T, N, V<:AbstractArray{T,1}}
		@assert length(v)==N
		new{T, N}(SVector{N, T}(v))
	end
end

merge(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}) where {T,K,L} = CamparyFloat{T, K+L}([x.val ; y.val]) 



end # module
