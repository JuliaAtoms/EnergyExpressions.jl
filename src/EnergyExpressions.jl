module EnergyExpressions

using AtomicLevels
using LinearAlgebra
using SparseArrays
using Combinatorics
using UnicodeFun
using Formatting
using ProgressMeter

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
    Base.:(:)(I::CartesianIndex{N}, J::CartesianIndex{N}) where N =
        CartesianIndices(map((i,j) -> i:j, Tuple(I), Tuple(J)))
end

include("conjugate_orbitals.jl")
include("slater_determinants.jl")
include("show_coefficients.jl")
include("nbody_operators.jl")
include("loop_macros.jl")
include("minors.jl")
include("nbody_matrix_elements.jl")
include("bit_configurations.jl")
include("common_operators.jl")
include("equations.jl")
include("variations.jl")
include("multi_configurational_equations.jl")
include("invariant_sets.jl")

end # module
