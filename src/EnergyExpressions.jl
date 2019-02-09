module EnergyExpressions

using AtomicLevels
using LinearAlgebra
using SparseArrays
using Combinatorics
using UnicodeFun

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
end

include("conjugate_orbitals.jl")
include("slater_determinants.jl")
include("nbody_operators.jl")
include("nbody_matrix_elements.jl")
include("common_operators.jl")
include("equations.jl")
include("variations.jl")
# include("invariant_sets.jl")

end # module
