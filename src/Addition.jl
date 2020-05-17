merge(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}) where {T,K,L} = CamparyFloat{K+L}([x.val ; y.val])

Base.:+(x::CamparyFloat{T,K}, y::CamparyFloat{T,L}, ::Val{M}) where {T,K,L, M} = Renorm(merge(x, y), Val(M)) 
Base.:+(x::CamparyFloat{T,K}, y::CamparyFloat{T,K}) where {T,K} = Base.:+(x, y, Val(K))