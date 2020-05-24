module CamparyBindings
using CuArrays


include("ErrorFree.jl")
include("Constructor.jl")
include("Addition.jl")
include("Renorm.jl")
include("PAddition.jl")

export CamparyFloat, FastTwoSum, FastTwoDiff, TwoSum, TwoProd, VecSum, VecSum!, FastVecSum, FastVecSum!

end # module
