using EnergyExpressions
import EnergyExpressions: @above_diagonal_loop, cofactor, nonzero_minors, find_diagonal_blocks,
    detaxis, powneg1, permutation_sign, indexsum,
    NBodyTerm, NBodyMatrixElement, OrbitalMatrixElement,
    contract,
    BitPattern, Orbitals, BitConfiguration,
    showbin, num_particles, orbital_numbers, get_orbitals,
    num_excitations, holes_mask, holes, particles_mask, particles,
    holes_particles_mask, holes_particles,
    relative_sign,
    orbital_overlap_matrix, non_zero_cofactors

using SparseArrays
using Combinatorics

using AtomicLevels
using LinearAlgebra
using Test
using DeepDiffs
using Suppressor

function strdisplay(o; kwargs...)
    buf = IOBuffer()
    io = IOContext(buf, kwargs...)
    show(io, MIME"text/plain"(), o)
    String(take!(buf))
end

rstripeachline(s) = join([rstrip(l) for l in split(s, '\n')], '\n')

function test_value_and_output(fun::Function, expected_value, expected_output)
    value = nothing
    out = @capture_out begin
        value = fun()
    end
    @test value == expected_value
    out = rstripeachline(out)
    expected_output = rstripeachline(expected_output)
    if out ≠ expected_output
        @error "Output does not match expected"
        println(deepdiff(expected_output, out))
    end
    @test out == expected_output
end

ome(a, op, b) = OrbitalMatrixElement(a, op, b)
ome(a::T, op, b::T) where {T<:Union{Symbol,<:Integer}} = ome([a], op, [b])
c(a, op, b) = ContractedOperator(Tuple(e for e in a), op, Tuple(e for e in b))
c(a::T, op, b::T) where {T<:Union{Symbol,<:Integer}} = c([a], op, [b])
eq(orb, op, args...) = NBodyEquation(orb, op, args...)
ceq(orb, op, args...) = eq(Conjugate(orb), op, args...)
ov(a,b) = OrbitalOverlap(a,b)
nbt(args...;f=1) = NBodyTerm([args...], f)
nbme(args...) = NBodyMatrixElement(NBodyTerm[args...])

function test_variations_equal(E, orbital, b)
    a = diff(E, orbital)
    a == b || @error "When varying, expected equality between" E a b a-b
    @test a == b
    b == a || @error "When varying, expected equality between" E b a b-a
    @test b == a
end

function test_equations_equal(a, b)
    a == b || @error "Expected equality between" a b a-b
    @test a == b
    b == a || @error "Expected equality between" b a b-a
    @test b == a
end

@testset "EnergyExpressions" verbose=true begin
    include("key_tracker.jl")
    include("locked_dict.jl")

    include("conjugate_orbitals.jl")
    include("nbody_matrix_elements.jl")
    include("bit_configurations.jl")
    include("loop_macros.jl")
    include("minors.jl")
    include("slater_determinants.jl")

    include("show_coefficients.jl")
    include("nbody_operators.jl")
    include("common_operators.jl")

    include("equations.jl")
    include("variations.jl")
end
