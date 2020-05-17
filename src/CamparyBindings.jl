module CamparyBindings
using StaticArrays

include("ErrorFree.jl")
include("Constructor.jl")
include("Addition.jl")
include("Renorm.jl")

export CamparyFloat, FastTwoSum, FastTwoDiff, TwoSum, TwoProd, VecSum, VecSum!, FastVecSum, FastVecSum!





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
