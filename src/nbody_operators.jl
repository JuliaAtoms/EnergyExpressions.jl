# * N-body operators

abstract type QuantumOperator end

"""
    NBodyOperator{N}

Abstract N-body operator coupling `N` bodies each between two Slater determinants
"""
abstract type NBodyOperator{N} <: QuantumOperator end

"""
    numbodies(::NBodyOperator{N})

Returns the number of bodies coupled by the N-body operator, i.e. `N`.
"""
numbodies(::NBodyOperator{N}) where N = N


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
    operators::Vector{<:Pair{<:NBodyOperator,<:Number}}
end

Base.zero(::LinearCombinationOperator) =
    LinearCombinationOperator(NBodyOperator[])
Base.iszero(op::LinearCombinationOperator) =
    isempty(op.operators)

"""
    numbodies(lco::LinearCombinationOperator)

Returns the maximum number of bodies coupled by any of the N-body
operators in the [`LinearCombinationOperator`](@ref).
"""
numbodies(lco::LinearCombinationOperator) =
    iszero(lco) ? 0 : maximum(numbodies.(first.(lco.operators)))

function Base.show(io::IO, op::LinearCombinationOperator)
    for (i,(o,n)) in enumerate(op.operators)
        i > 1 && write(io, " ")
        showcoeff(io, n, i > 1)
        write(io, string(o))
    end
end

Base.:(+)(a::NBodyOperator, b::NBodyOperator) =
    LinearCombinationOperator([a=>1,b=>1])
Base.:(-)(a::NBodyOperator, b::NBodyOperator) =
    LinearCombinationOperator([a=>1,b=>-1])

Base.:(+)(a::LinearCombinationOperator, b::NBodyOperator) =
    LinearCombinationOperator(vcat(a.operators, b=>1))
Base.:(-)(a::LinearCombinationOperator, b::NBodyOperator) =
    LinearCombinationOperator(vcat(a.operators, b=>-1))

Base.:(+)(a::NBodyOperator, b::LinearCombinationOperator) =
    LinearCombinationOperator(vcat(a=>1, b.operators))
Base.:(-)(a::NBodyOperator, b::LinearCombinationOperator) =
    LinearCombinationOperator(vcat(a=>1, Pair{<:NBodyOperator,<:Number}[o[1]=>-o[2]
                                                                       for o in b.operators]))

Base.:(+)(a::LinearCombinationOperator, b::LinearCombinationOperator) =
    LinearCombinationOperator(vcat(a.operators, b.operators))
Base.:(-)(a::LinearCombinationOperator, b::LinearCombinationOperator) =
    LinearCombinationOperator(vcat(a.operators, Pair{<:NBodyOperator,<:Number}[o[1]=>-o[2]
                                                                               for o in b.operators]))

Base.:(*)(a::Number, b::NBodyOperator) =
    LinearCombinationOperator([b=>a])
Base.:(*)(a::NBodyOperator, b::Number) =
    LinearCombinationOperator([a=>b])

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

Base.:(==)(a::ContractedOperator, b::ContractedOperator) =
    a.a == b.a && a.o == b.o && a.b == b.b

"""
    in(orbital, co::ContractedOperator)

Test if `orbital` is among the right set of orbitals of the
[`ContractedOperator`](@ref) `co`. Useful to test if `co` is an
integral operator with respect to `orbital`.
"""
Base.in(orbital::O, co::ContractedOperator) where O = orbital ∈ co.b
"""
    in(corbital::Conjugate, co::ContractedOperator)

Test if `corbital` is among the left set of orbitals of the
[`ContractedOperator`](@ref) `co`. Useful to test if `co` is an
integral operator with respect to `corbital`.
"""
Base.in(corbital::Conjugate{O}, co::ContractedOperator) where O = corbital.orbital ∈ co.a

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
complement(N::Integer, i::Integer...) = setdiff(1:N, i)
complement(N::Integer, i::AbstractVector{<:Integer}) = setdiff(1:N, i)

function Base.show(io::IO, o::ContractedOperator{N}) where N
    write(io, "⟨", join(string.(o.a), " "))
    write(io, "|")
    show(io, o.o)
    write(io, "|")
    write(io, join(string.(o.b), " "), "⟩")
end

export QuantumOperator, NBodyOperator, numbodies,
    ZeroBodyOperator, OneBodyOperator, TwoBodyOperator,
    IdentityOperator, ContractedOperator
