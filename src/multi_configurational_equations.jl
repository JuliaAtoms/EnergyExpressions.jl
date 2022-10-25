"""
    MCTerm(i, j, coeff, operator, source_orbital, integrals=[])

Represents one term in the multi-configurational expansion. `i` and
`j` are indices in the mixing-coefficient vector c (which is subject
to optimization, and thus has to be referred to), `coeff` is an
additional coefficient, and `integrals` is a list of indices into the
vector of common integrals, the values of which should be multiplied
to form the overall coefficient.
"""
struct MCTerm{T,QO,O}
    i::Int
    j::Int
    coeff::T
    operator::QO
    source_orbital::O
    integrals::Vector{Int}
end

Base.:(==)(a::MCTerm, b::MCTerm) =
    a.i == b.i && a.j == b.j &&
    a.coeff == b.coeff && a.operator == b.operator &&
    a.source_orbital == b.source_orbital &&
    sort(a.integrals) == sort(b.integrals)

"""
    OrbitalEquation(orbital, equation,
                    one_body, direct_terms, exchange_terms, source_terms)

Represents the integro-differential equation for `orbital`, expressed
as a linear combination of the different terms, with pointers to the
list of common integrals that is stored by the encompassing
[`MCEquationSystem`](@ref) object.
"""
struct OrbitalEquation{O,Equation}
    orbital::O
    equation::Equation
    terms::Vector{Pair{Int,Vector{MCTerm}}}
end

Base.:(==)(a::OrbitalEquation, b::OrbitalEquation) =
    a.orbital == b.orbital &&
    a.equation == b.equation &&
    a.terms == b.terms


function Base.show(io::IO, oeq::OrbitalEquation)
    write(io, "OrbitalEquation($(oeq.orbital)): ")
    show(io, oeq.equation)
    write(io, "\n")
end

Base.iszero(oeq::OrbitalEquation) = iszero(oeq.equation)

"""
    MCEquationSystem(equations, integrals)

Represents a coupled system of integro-differential `equations`,
resulting from the variation of a multi-configurational
[`EnergyExpression`](@ref), with respect to all constituent
orbitals. All `integrals` that are in common between the `equations`
need only be computed once per iteration, for efficiency.
"""
struct MCEquationSystem
    equations
    integrals
end

function Base.show(io::IO, mceqs::MCEquationSystem)
    neq = length(mceqs.equations)
    nnzeq = count(!iszero, mceqs.equations)
    nint = length(mceqs.integrals)

    write(io, "$(neq)-equation MCEquationSystem, $(nnzeq) non-zero equations, $(nint) common integrals")
end

function Base.show(io::IO, mime::MIME"text/plain", mceqs::MCEquationSystem)
    show(io, mceqs)
    println(io)
    nd = length(digits(length(mceqs.equations)))
    for (i,eq) in enumerate(mceqs.equations)
        iszero(eq) && continue
        printfmt(io, "- {1:>$(nd)d}: ", i)
        show(io, mime, eq)
    end
end

"""
    pushifmissing!(vector, element)

Push `element` to the end of `vector`, if not already present. Returns
the index of `element` in `vector`.
"""
function pushifmissing!(v::Vector, element)
    i = findfirst(isequal(element), v)
    isnothing(i) || return i
    push!(v, element)
    length(v)
end

"""
    orbital_equation(E::EnergyExpression, orbital, integrals::Vector)

Generate the [`OrbitalEquation`](@ref) governing `orbital` by varying
the [`EnergyExpression`](@ref) `E`, and storing common expressions in
`integrals`.
"""
function orbital_equation(E::EM, orbital::O, integrals::KeyTracker) where {EM<:EnergyExpression,O}
    equation = diff(E, orbital)

    # Complementary orbital
    comp_orbital = orbital isa Conjugate ? orbital.orbital : Conjugate(orbital)

    terms = Dict{Int,Vector{MCTerm}}()

    # Invert matrix coordinate -> equation mapping, i.e. gather all
    # coordinates for which a specific NBodyEquation appears.

    for (i,j,eq) in zip(findnz(equation)...)
        for subeq = eq.equations
            operator = subeq.operator
            source_orbital = subeq.orbital
            coeff = MCTerm(i,j,subeq.factor.coeff,
                           operator, source_orbital,
                           map(factor -> get!(integrals, factor),
                               subeq.factor.factors) |> Vector{Int})

            # If the operator isa ContractedOperator, an integral has
            # to be performed, which can potentially be reused among
            # common terms. This is only possible, however, if the
            # orbital under consideration is not inside the integrand;
            # if it is, we are dealing with an integral operator,
            # which has to be reevaluated each time it is applied to
            # an orbital.
            integral = (operator isa ContractedOperator && comp_orbital ∉ operator) ?
                get!(integrals, operator) : 0
            terms[integral] = push!(get(terms, integral, MCTerm[]), coeff)
        end
    end

    OrbitalEquation(comp_orbital, equation, collect(pairs(terms)))
end

"""
    diff(energy_expression, orbitals)

Derive the integro-differential equations for all `orbitals`, from
`energy_expression`. Returns a [`MCEquationSystem`](@ref), that
gathers information on which integrals are common to all equations,
for efficient equation solving.
"""
function Base.diff(E::EM, orbitals::VO; verbosity=0) where {EM<:EnergyExpression, O,VO<:AbstractVector{O}}
    # Vector of common integrals
    integrals = KeyTracker(LockedDict{Any,Int}())

    norb = length(orbitals)
    p = if verbosity > 0
        @info "Deriving equations for $(norb) orbitals"
        Progress(norb)
    end

    equations = [OrbitalEquation[] for i in 1:Threads.nthreads()]
    for i = 1:length(orbitals)
        tid = Threads.threadid()
        orbital = orbitals[i]
        eq = orbital_equation(E, orbital, integrals)
        isnothing(p) || ProgressMeter.next!(p)
        push!(equations[tid], eq)
    end
    equations = reduce(vcat, equations)

    MCEquationSystem(equations, keys(integrals))
end

"""
    diff(fun!, energy_expression, orbitals)

Derive the integro-differential equations for all `orbitals`, from
`energy_expression`; after each orbital equation has been generated
`fun!` is applied to `energy_expression` with the current orbital as
argument, which allows gradual modification of the energy
expression. Returns a [`MCEquationSystem`](@ref), that gathers
information on which integrals are common to all equations, for
efficient equation solving.

"""
function Base.diff(fun!::Function, E::EM, orbitals::VO; verbosity=0) where {EM<:EnergyExpression, O,VO<:AbstractVector{O}}
    # Vector of common integrals
    integrals = KeyTracker{Any}()

    E = deepcopy(E)

    norb = length(orbitals)
    p = if verbosity > 0
        @info "Deriving equations for $(norb) orbitals"
        Progress(norb)
    end
    equations = map(orbitals) do orbital
        eq = orbital_equation(E, orbital, integrals)
        fun!(E, orbital)
        isnothing(p) || ProgressMeter.next!(p)
        eq
    end

    MCEquationSystem(equations, keys(integrals))
end
