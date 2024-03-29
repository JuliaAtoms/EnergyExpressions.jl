# * N-body operators

abstract type QuantumOperator end

"""
    ishermitian(op::QuantumOperator)

By default, all `QuantumOperator`s are Hermitian; this can be
overridden for subtypes to explicitly declare an operator
non-Hermitian.
"""
LinearAlgebra.ishermitian(::QuantumOperator) = true

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

LinearAlgebra.adjoint(I::IdentityOperator) = I

# * Adjoint operator

struct AdjointOperator{N,O<:NBodyOperator{N}} <: NBodyOperator{N}
    o::O
end

Base.parent(ao::AdjointOperator) = ao.o

Base.hash(ao::AdjointOperator, h::UInt) = hash(parent(ao), hash(0x123, h))
numbodies(ao::AdjointOperator) = numbodies(parent(ao))

LinearAlgebra.adjoint(o::QuantumOperator) = AdjointOperator(o)
LinearAlgebra.adjoint(ao::AdjointOperator) = parent(ao)

function Base.show(io::IO, ao::AdjointOperator)
    show(io, parent(ao))
    write(io, "†")
end

# * Linear combination

"""
    LinearCombinationOperator(operators)

Represents a linear combination of [`NBodyOperator`](@ref)s.
"""
struct LinearCombinationOperator <: QuantumOperator
    operators::Vector{<:Pair{<:NBodyOperator,<:Number}}
end

Base.zero(::LinearCombinationOperator) =
    LinearCombinationOperator(Pair{<:NBodyOperator,<:Number}[])
Base.zero(::Type{LinearCombinationOperator}) =
    LinearCombinationOperator(Pair{<:NBodyOperator,<:Number}[])

Base.iszero(op::LinearCombinationOperator) =
    isempty(op.operators)

LinearAlgebra.adjoint(op::LinearCombinationOperator) =
    sum(adjoint(o)*conj(c) for (o,c) in op.operators)

Base.:(==)(a::LinearCombinationOperator,b::LinearCombinationOperator) =
    a.operators == b.operators

function Base.hash(lco::LinearCombinationOperator, h::UInt)
    for (op,coeff) in lco.operators
        h = hash(op, hash(coeff, h))
    end
    h
end

"""
    numbodies(lco::LinearCombinationOperator)

Returns the maximum number of bodies coupled by any of the N-body
operators in the [`LinearCombinationOperator`](@ref).
"""
numbodies(lco::LinearCombinationOperator) =
    iszero(lco) ? 0 : maximum(numbodies.(first.(lco.operators)))

function Base.show(io::IO, op::LinearCombinationOperator)
    if iszero(op)
        write(io, "0")
        return
    end
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

Base.:(-)(op::QuantumOperator) = -1op

Base.:(*)(a::Number, b::LinearCombinationOperator) =
    LinearCombinationOperator(Pair{<:NBodyOperator,<:Number}[op=>a*c for (op,c) in b.operators])
Base.:(*)(a::LinearCombinationOperator, b::Number) =
    LinearCombinationOperator(Pair{<:NBodyOperator,<:Number}[op=>b*c for (op,c) in a.operators])

Base.filter(fun::Function, op::LinearCombinationOperator) =
    LinearCombinationOperator(filter(oc -> fun(oc[1]), op.operators))

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
function contract end

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

export QuantumOperator, NBodyOperator, AdjointOperator, numbodies,
    ZeroBodyOperator, OneBodyOperator, TwoBodyOperator,
    IdentityOperator, ContractedOperator
