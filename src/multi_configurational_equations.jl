"""
    MCCoeff(i, j, coeff, integrals=[])

Represents the coefficient of one term in the multi-configurational
expansion. `i` and `j` are indices in the mixing-coefficient vector c
(which is subject to optimization, and thus has to be referred to),
`coeff` is an additional coefficient, and `integrals` is a list of
indices into the vector of common integrals, the values of which
should be multiplied to form the overall coefficient.
"""
struct MCCoeff{T}
    i::Int
    j::Int
    coeff::T
    integrals::Vector{Int}
end

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
    one_body::Vector{MCCoeff}
    direct_terms::Vector{Pair{Int,Vector{MCCoeff}}}
    exchange_terms::Vector{Tuple{Any,MCCoeff,Any}}
    source_terms::Vector{Pair{Int,Vector{Tuple{MCCoeff,Any}}}}
end

function Base.show(io::IO, oeq::OrbitalEquation)
    write(io, "OrbitalEquation($(oeq.orbital)): ")
    show(io, oeq.equation)
    write(io, "\n")
end

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
function orbital_equation(E::EM, orbital::O, integrals::Vector) where {EM<:EnergyExpression,O}
    equation = diff(E, orbital)

    # Complementary orbital
    comp_orbital = orbital isa Conjugate ? orbital.orbital : Conjugate(orbital)

    # For each orbital, there are a few different kinds of terms
    # possible:
    # - Acting on the orbital in question:
    #   - One-body Hamiltonian
    #   - Direct potentials
    #   - Exchange potentials
    # - Source terms, resulting from the configuration interaction
    one_body = MCCoeff[]
    direct_terms = Dict{Int,Vector{MCCoeff}}()
    exchange_terms = Tuple{Any,MCCoeff,Any}[]
    source_terms = Dict{Int,Vector{Tuple{MCCoeff,Any}}}()

    # Invert matrix coordinate -> equation mapping, i.e. gather all
    # coordinates for which a specific NBodyEquation appears.

    for (i,j,eq) in zip(findnz(equation)...)
        for subeq = eq.equations
            operator = subeq.operator
            source_orbital = subeq.orbital
            coeff = MCCoeff(i,j,subeq.factor.coeff,
                            map(factor -> pushifmissing!(integrals, factor),
                                subeq.factor.factors) |> Vector{Int})

            if operator isa FieldFreeOneBodyHamiltonian
                source_orbital == comp_orbital ||
                    throw(ArgumentError("Field-free one-body Hamiltonian not pertaining to $(orbital)"))
                push!(one_body, coeff)
            elseif operator isa ContractedOperator && (comp_orbital ∈ operator.b ||
                                                      source_orbital == comp_orbital)
                if comp_orbital ∈ operator.b
                    push!(exchange_terms, (operator,coeff,source_orbital))
                else
                    integral = pushifmissing!(integrals, operator)
                    direct_terms[integral] =
                        push!(get(direct_terms, integral, MCCoeff[]),
                              coeff)
                end
            else
                integral = pushifmissing!(integrals, operator)
                source_terms[integral] =
                    push!(get(source_terms, integral, Tuple{MCCoeff,Any}[]),
                          (coeff, source_orbital))
            end
        end
    end

    OrbitalEquation(comp_orbital, equation, one_body,
                    collect(pairs(direct_terms)), exchange_terms,
                    collect(pairs(source_terms)))
end

"""
    diff(energy_expression, orbitals)

Derive the integro-differential equations for all `orbitals`, from
`energy_expression`. Returns a [`MCEquationSystem`](@ref), that
gathers information on which integrals are common to all equations,
for efficient equation solving.
"""
function Base.diff(E::EM, orbitals::VO) where {EM<:EnergyExpression, O,VO<:AbstractVector{O}}
    # Vector of common integrals
    integrals = []

    equations = map(orbitals) do orbital
        orbital_equation(E, orbital, integrals)
    end

    MCEquationSystem(equations, integrals)
end
