# * N-body operators

abstract type QuantumOperator end

"""
    NBodyOperator{N}

Abstract N-body operator coupling `N` bodies each between two Slater determinants
"""
abstract type NBodyOperator{N} <: QuantumOperator end


# # Examples

# ```jldoctest
# julia> ZeroBodyOperator
# Ω₀

# julia> OneBodyOperator
# Ω₁

# julia> TwoBodyOperator
# Ω₂

# julia> NBodyOperator{6}
# Ω₆
# ```

# Base.show(io::IO, ::Type{NBodyOperator{N}}) where N =
#     write(io, "Ω", to_subscript(N))

const ZeroBodyOperator = NBodyOperator{0}
const OneBodyOperator = NBodyOperator{1}
const TwoBodyOperator = NBodyOperator{2}

"""
    IdentityOperator{N}

The N-body identity operator. Leaves the orbital(s) acted upon
unchanged.
"""
struct IdentityOperator{N} <: NBodyOperator{N} end

Base.show(io::IO, ::IdentityOperator{N}) where N =
    write(io, to_boldface("I"), to_subscript(N))

# * Linear combination

"""
    LinearCombinationOperator(operators)

Represents a linear combination of [`NBodyOperator`](@ref)s.
"""
struct LinearCombinationOperator <: QuantumOperator
    operators::Vector{<:NBodyOperator}
end

Base.zero(::LinearCombinationOperator) =
    LinearCombinationOperator(NBodyOperator[])
Base.iszero(op::LinearCombinationOperator) =
    isempty(op.operators)

Base.show(io::IO, op::LinearCombinationOperator) =
    write(io, join(string.(op.operators), " + "))

Base.:(+)(a::NBodyOperator, b::NBodyOperator) =
    LinearCombinationOperator(NBodyOperator[a,b])

Base.:(+)(a::LinearCombinationOperator, b::NBodyOperator) =
    LinearCombinationOperator(vcat(a.operators, b))
Base.:(+)(a::NBodyOperator, b::LinearCombinationOperator) =
    LinearCombinationOperator(vcat(a, b.operators))

# * Contracted operators

"""
    ContractedOperator(a, o, b)

An [`NBodyOperator`](@ref) representing the contraction of the
operator `o` over the orbital sets `a` and `b`. The lengths of `a` and
`b` have to equal, and they cannot exceed the dimension of `o`.
"""
struct ContractedOperator{N,P,Q,A,O<:NBodyOperator{P},B} <: NBodyOperator{N}
    a::NTuple{Q,A}
    o::O
    b::NTuple{Q,B}
    function ContractedOperator(a::NTuple{Q,A}, o::O, b::NTuple{Q,B}) where {P,Q,A,O<:NBodyOperator{P},B}
        P ≥ Q || throw(ArgumentError("Cannot contract $(P)-body operator over $(Q) coordinates"))
        N = P - Q
        new{N,P,Q,A,O,B}(a, o, b)
    end
end

"""
    contract(orbital_matrix_element, i...)

Contract the `orbital_matrix_element` over all coordinates `i...`.
"""
contract(ome, i::Integer...) =
    throw(ArgumentError("`contract` not implemented for $(ome)"))

"""
    complement(N, i...)

Generate the complement to `i...` in the set `1:N`. Useful for
contracting [`OrbitalMatrixElement`](@ref)s over all coordinates
_except_ `i...`.
"""
function complement(N::Integer, i::Integer...)
    i = sort([i...])
    filter(!isempty,
           vcat(1:i[1]-1,
                [(i[j-1]+1):(i[j]-1) for j in 2:length(i)]...,
                i[end]+1:N))
end

function Base.show(io::IO, o::ContractedOperator{N}) where N
    write(io, "⟨", join(string.(o.a), " "))
    write(io, "|")
    show(io, o.o)
    write(io, "|")
    write(io, join(string.(o.b), " "), "⟩")
end

export NBodyOperator, ZeroBodyOperator, OneBodyOperator, TwoBodyOperator, ContractedOperator
