# * Slater determinants

"""
    SlaterDeterminant(orbitals::Vector{O})

Constructs a Slater determinant from a set of spin-orbitals.

# Examples

```jldoctest
julia> SlaterDeterminant([:a, :b])
a(1)b(2) - a(2)b(1)

julia> SlaterDeterminant([:a, :b, :c])
a(1)b(2)c(3) - a(1)b(3)c(2) - a(2)b(1)c(3) + a(2)b(3)c(1) + a(3)b(1)c(2) - a(3)b(2)c(1)
```
"""
struct SlaterDeterminant{O}
    orbitals::Vector{O}
end

"""
    SlaterDeterminant(configuration::Configuration{<:SpinOrbital})

Constructs a Slater determinant from the spin-orbitals of the
spin-configuration `configuration`.

# Examples

```jldoctest
julia> SlaterDeterminant(spin_configurations(c"1s2")[1])
1s₀α(1)1s₀β(2) - 1s₀α(2)1s₀β(1)
```
"""
SlaterDeterminant(c::Configuration{<:SpinOrbital}) =
    SlaterDeterminant(c.orbitals)

"""
    length(slater_determinant)

Return the number of spin-orbitals in the Slater determinant.
"""
Base.length(sd::SlaterDeterminant) = length(sd.orbitals)

function Base.show(io::IO, sd::SlaterDeterminant)
    n = length(sd)
    for (j,p) in enumerate(permutations(1:n))
        s = Combinatorics.parity(p)
        (j > 1 || s == 1) && write(io, " ", s == 0 ? "+" : "-", " ")
        for (i,c) in enumerate(p)
            write(io, "$(sd.orbitals[i])($(c))")
        end
    end
end

"""
    AdjointSlaterDeterminant(slater_determinant)

Representation of the Hermitian conjugate (dual vector) of a Slater
determinant. Constructed via the usual [`adjoint`](@ref) operator.
"""
struct AdjointSlaterDeterminant{O}
    sd::SlaterDeterminant{O}
end

"""
    adjoint(slater_determinant)

Construct the adjoint of `slater_determinant`

# Examples

```jldoctest
julia> SlaterDeterminant([:a, :b])'
[a(1)b(2) - a(2)b(1)]†
```
"""
Base.adjoint(sd::SlaterDeterminant) =
    AdjointSlaterDeterminant(sd)

"""
    length(adjoint_slater_determinant)

Return the number of spin-orbitals in the adjoint Slater determinant.
"""
Base.length(asd::AdjointSlaterDeterminant) = length(asd.sd)

function Base.show(io::IO, asd::AdjointSlaterDeterminant)
    write(io, "[")
    show(io, asd.sd)
    write(io, "]†")
end

# * N-body operators

"""
    NBodyOperator{N}

N-body operator coupling `N` bodies each between two Slater determinants

# Examples

```jldoctest
julia> ZeroBodyOperator
Ω₀

julia> OneBodyOperator
Ω₁

julia> TwoBodyOperator
Ω₂

julia> NBodyOperator{6}
Ω₆
```
"""
struct NBodyOperator{N} end

Base.show(io::IO, ::Type{NBodyOperator{N}}) where N =
    write(io, "Ω", to_subscript(N))

const ZeroBodyOperator = NBodyOperator{0}
const OneBodyOperator = NBodyOperator{1}
const TwoBodyOperator = NBodyOperator{2}

# ** Matrix elements

"""
    NBodyTermFactor

Abstract type for a factor in a term in a N-body matrix element expansion
"""
abstract type NBodyTermFactor end

"""
    OrbitalOverlap(a,b)

Represents the overlap between the orbitals `a` and `b` in a N-body matrix element expansion.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalOverlap(:a,:b)
⟨a|b⟩
```
"""
struct OrbitalOverlap{A,B} <: NBodyTermFactor
    a::A
    b::B
end

Base.show(io::IO, o::OrbitalOverlap) =
    write(io, "⟨$(o.a)|$(o.b)⟩")


"""
    OrbitalMatrixElement(a,o,b)

Represents the N-body matrix element between the sets of orbitals `a` and `b`.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), TwoBodyOperator, (:c,:d))
⟨a b|Ω₂|c d⟩
```
"""
struct OrbitalMatrixElement{N,A,O,B} <: NBodyTermFactor
    a::NTuple{N,A}
    o::O
    b::NTuple{N,B}    
end

function Base.show(io::IO, ome::OrbitalMatrixElement{N}) where N
    write(io, "⟨", join(string.(ome.a), " "))
    N > 0 && write(io, "|")
    show(io, ome.o)
    N > 0 && write(io, "|")
    write(io, join(string.(ome.b), " "), "⟩")
end

"""
    NBodyTerm(factors, sign)

Structure representing one term in the expansion of a N-body matrix element.
"""
struct NBodyTerm
    factors::Vector{NBodyTermFactor}
    sign
end

function Base.show(io::IO, term::NBodyTerm)
    write(io, term.sign == -1 ? "- " : "+ ")
    write(io, join(string.(term.factors), ""))
end

"""
    NBodyMatrixElement(terms)

Structure representing the expansion of a N-body matrix element.
"""
struct NBodyMatrixElement
    terms::Vector{NBodyTerm}
end

Base.show(io::IO, nbme::NBodyMatrixElement) =
    write(io, join(string.(nbme.terms), " "))

function Base.:(*)(a::AdjointSlaterDeterminant, ::Type{ZeroBodyOperator}, b::SlaterDeterminant)
    length(a) == length(b) ||
        throw(DimensionMismatch("Slater determinants of different dimensions"))

    aorbitals = a.sd.orbitals
    borbitals = b.orbitals
    
    ps = permutations(1:length(a))
    
    terms = map(ps) do pa
        sa = Combinatorics.parity(pa) == 1 ? -1 : 1
        map(ps) do pb
            sb = Combinatorics.parity(pb) == 1 ? -1 : 1
            factors = map(zip(pa,pb)) do (ai,bi)
                OrbitalOverlap(aorbitals[ai],borbitals[bi])
            end
            NBodyTerm(factors, sa*sb)
        end |> t -> vcat(t...)
    end |> t -> vcat(t...)

    NBodyMatrixElement(terms)
end

function Base.:(*)(a::AdjointSlaterDeterminant, ::Type{OneBodyOperator}, b::SlaterDeterminant)
    length(a) == length(b) ||
        throw(DimensionMismatch("Slater determinants of different dimensions"))

    aorbitals = a.sd.orbitals
    borbitals = b.orbitals

    n = length(a)
    ps = permutations(1:n)
    
    terms = map(ps) do pa
        sa = Combinatorics.parity(pa) == 1 ? -1 : 1
        map(ps) do pb
            sb = Combinatorics.parity(pb) == 1 ? -1 : 1
            factors = map(zip(pa,pb)) do (ai,bi)
                OrbitalOverlap(aorbitals[ai],borbitals[bi])
            end
            map(1:n) do i
                op = OrbitalMatrixElement((aorbitals[pa[i]],),
                                          OneBodyOperator,
                                          (borbitals[pb[i]],))
                NBodyTerm(vcat(factors[1:i-1], op, factors[i+1:end]), sa*sb)
            end
        end |> t -> vcat(t...)
    end |> t -> vcat(t...)

    NBodyMatrixElement(terms)
end

struct NBodyProxy{ASD<:AdjointSlaterDeterminant,O}
    asd::ASD
    o::O
end

Base.:(*)(asd::AdjointSlaterDeterminant, ::Type{NBodyOperator{N}}) where {N} =
    NBodyProxy(asd, NBodyOperator{N})

Base.:(*)(nbp::NBodyProxy, sd::SlaterDeterminant) =
    *(nbp.asd, nbp.o, sd)


export SlaterDeterminant, NBodyOperator, ZeroBodyOperator, OneBodyOperator, TwoBodyOperator
