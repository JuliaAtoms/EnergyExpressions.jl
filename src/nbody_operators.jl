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

export NBodyOperator, ZeroBodyOperator, OneBodyOperator, TwoBodyOperator
