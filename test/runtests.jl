using EnergyExpressions
using AtomicLevels
using LinearAlgebra
using Test

function strdisplay(o)
    buf = IOBuffer()
    show(buf, MIME"text/plain"(), o)
    String(take!(buf))
end

include("key_tracker.jl")
include("locked_dict.jl")

include("conjugate_orbitals.jl")
include("bit_configurations.jl")
include("loop_macros.jl")
include("minors.jl")
include("slater_determinants.jl")

include("show_coefficients.jl")
include("common_operators.jl")

include("variations.jl")
