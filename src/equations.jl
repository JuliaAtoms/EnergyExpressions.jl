"""
    NBodyEquation{N,O}(orbital, operator::NBodyOperator[, factor::NBodyTerm])

Equation for an `orbital`, acted upon by an operator, which may be a
single-particle operator, or an N-body operator, contracted over all
coordinates but one, and optionally multiplied by an
[`NBodyTerm`](@ref), corresponding to overlaps/matrix elements between
other orbitals.
"""
struct NBodyEquation{QO<:OneBodyOperator,O}
    orbital::O
    operator::QO
    factor::NBodyTerm
    NBodyEquation(orbital::O, operator::QO,
                  factor::NBodyTerm=one(NBodyTerm)) where {QO<:OneBodyOperator,O} =
        new{QO,O}(orbital, operator, factor)
end

Base.iszero(::NBodyEquation) = false

Base.:(*)(nbe::NBodyEquation, factor::NBodyTerm) =
    NBodyEquation(nbe.orbital, nbe.operator, nbe.factor*factor)

function Base.show(io::IO, eq::NBodyEquation)
    show(io, eq.factor, show_sign=true)
    eq.orbital isa Conjugate && write(io, "⟨$(conj(eq.orbital))|")
    show(io, eq.operator)
    !(eq.orbital isa Conjugate) && write(io, "|$(eq.orbital)⟩")
end

"""
    LinearCombinationEquation(equations)

A type representing a linear combination of
[`NBodyEquation`](@ref)s. Typically arises when varying a multi-term
energy expression.
"""
struct LinearCombinationEquation
    equations::Vector{<:NBodyEquation}
end

Base.zero(::Type{LinearCombinationEquation}) = LinearCombinationEquation(NBodyEquation[])
Base.zero(::LinearCombinationEquation) = zero(LinearCombinationEquation)
Base.iszero(eq::LinearCombinationEquation) =
    isempty(eq.equations)

function Base.show(io::IO, eq::LinearCombinationEquation)
    if iszero(eq)
        write(io, "0")
        return
    end
    noprintterms = 2

    terms = if get(io, :limit, false) && length(eq.equations) > 2noprintterms
        vcat(eq.equations[1:noprintterms], "…", eq.equations[end-(noprintterms-1):end])
    else
        eq.equations
    end
    write(io, join(string.(terms), " "))
end

export NBodyEquation, LinearCombinationEquation
