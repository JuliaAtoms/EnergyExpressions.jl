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

Base.iszero(nbe::NBodyEquation) = iszero(nbe.factor)
Base.zero(::Type{NBodyEquation}) = NBodyEquation(Symbol(:vac), IdentityOperator{1}(), zero(NBodyTerm))

Base.:(*)(nbe::NBodyEquation, factor) =
    NBodyEquation(nbe.orbital, nbe.operator, nbe.factor*factor)

Base.:(*)(factor, nbe::NBodyEquation) = nbe*factor

Base.:(-)(nbe::NBodyEquation) = -1*nbe

function Base.show(io::IO, eq::NBodyEquation)
    show(io, eq.factor, show_sign=true)
    eq.orbital isa Conjugate && write(io, "⟨$(conj(eq.orbital))|")
    show(io, eq.operator)
    !(eq.orbital isa Conjugate) && write(io, "|$(eq.orbital)⟩")
end

function Base.:(+)(a::NBodyEquation, b::NBodyEquation)
    factor = a.factor + b.factor
    if a.orbital == b.orbital && a.operator == b.operator && factor isa NBodyTerm
        NBodyEquation(a.orbital, a.operator, factor)
    else
        LinearCombinationEquation([a]) + b
    end
end

Base.:(==)(a::NBodyEquation, b::NBodyEquation) =
    a.orbital == b.orbital && a.operator == b.operator && a.factor == b.factor

function add_equations!(equations::AbstractVector{<:NBodyEquation}, eq::NBodyEquation)
    i = findfirst(a -> a.orbital == eq.orbital && a.operator == eq.operator,
                  equations)
    if isnothing(i)
        push!(equations, eq)
    else
        factor = eq.factor+equations[i].factor
        if iszero(factor)
            deleteat!(equations, i)
        elseif factor isa NBodyTerm
            # If the resulting multiplicative factor is something
            # simple, like a plain number, we can combine the terms
            # into one.
            equations[i] = NBodyEquation(eq.orbital, eq.operator, factor)
        else
            # However, if we end with something complicated —
            # e.g. variation of R⁰(a,a;a,b)⟨a|a⟩ + R⁰(a,a;a,b) with
            # respect to ⟨a| would yield five terms
            #
            # + ⟨a|a⟩r⁻¹×Y⁰(a,b)|a⟩ + ⟨a|a⟩r⁻¹×Y⁰(a,a)|b⟩ + R⁰(a,a;a,b)𝐈₁|a⟩ + 1r⁻¹×Y⁰(a,b)|a⟩ + 1r⁻¹×Y⁰(a,a)|b⟩,
            #
            # two of which we could group as
            #
            # + 1r⁻¹×Y⁰(a,b)|a⟩ (1 + ⟨a|a⟩)
            #
            # — we simply treat them as separate terms. We do this
            # because it's easier to implement, but also because we
            # probably want to handle the terms separately (or do
            # we?). In any case, this should be pretty rare.
            push!(equations, eq)
        end
    end
    equations
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
Base.iszero(eq::LinearCombinationEquation) = isempty(eq.equations) || all(iszero, eq.equations)

Base.convert(::Type{LinearCombinationEquation}, eq::NBodyEquation) =
    LinearCombinationEquation([eq])

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

function Base.:(*)(eq::LinearCombinationEquation, factor)
    equations = [NBodyEquation(nbe.orbital, nbe.operator, nbe.factor*factor)
                 for nbe in eq.equations]
    LinearCombinationEquation(equations)
end

Base.:(*)(factor, eq::LinearCombinationEquation) = eq*factor
Base.:(-)(eq::LinearCombinationEquation) = -1*eq

function Base.:(+)(a::LinearCombinationEquation, b::NBodyEquation)
    T = promote_type(eltype(a.equations), typeof(b))
    LinearCombinationEquation(add_equations!(Vector{T}(copy(a.equations)), b))
end

Base.:(+)(a::NBodyEquation, b::LinearCombinationEquation) = b + a

function add_equations!(equations::AbstractVector{<:NBodyEquation}, eqs::LinearCombinationEquation)
    foreach(Base.Fix1(add_equations!, equations), eqs.equations)
    equations
end

function Base.:(+)(a::LinearCombinationEquation, b::LinearCombinationEquation)
    T = promote_type(eltype(a.equations), eltype(b.equations))
    LinearCombinationEquation(add_equations!(Vector{T}(copy(a.equations)), b))
end

Base.:(==)(a::LinearCombinationEquation, b::LinearCombinationEquation) =
    compare_vectors(a.equations, b.equations)

function Base.:(==)(a::LinearCombinationEquation, b::NBodyEquation)
    length(a.equations) == 0 && iszero(b) && return true
    length(a.equations) == 1 || return false
    a.equations[1] == b
end

Base.:(==)(a::NBodyEquation, b::LinearCombinationEquation) = (b == a)

function Base.:(-)(a::LinearCombinationEquation, b::LinearCombinationEquation)
    iszero(b) && return a
    iszero(a) && return -b
    eqs = Vector{NBodyEquation}(copy(a.equations))
    for e in b.equations
        i = findfirst(ae -> ae.orbital == e.orbital && ae.operator == e.operator, eqs)
        if isnothing(i)
            push!(eqs, -e)
        else
            factor = eqs[i].factor - e.factor
            if iszero(factor)
                deleteat!(eqs, i)
            elseif factor isa NBodyTerm
                eqs[i] = NBodyEquation(e.orbital, e.operator, factor)
            else
                push!(eqs, -e)
            end
        end
    end
    LinearCombinationEquation(eqs)
end

Base.:(-)(a::LinearCombinationEquation, b::NBodyEquation) = a - LinearCombinationEquation([b])
Base.:(-)(a::NBodyEquation, b::LinearCombinationEquation) = LinearCombinationEquation([a]) - b

export NBodyEquation, LinearCombinationEquation
