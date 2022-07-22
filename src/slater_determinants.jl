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
AtomicLevels.num_electrons(sd::SlaterDeterminant) = length(sd)
AtomicLevels.orbitals(sd::SlaterDeterminant) = sd.orbitals

Base.show(io::IO, sd::SlaterDeterminant) =
    write(io, "|",join(sd.orbitals, " "), "|")

function Base.show(io::IO, ::MIME"text/plain", sd::SlaterDeterminant)
    noprintterms = 5
    n = length(sd)
    if n > 3
        show(io, sd)
        return
    end
    numperm = factorial(n)
    for (j,p) in enumerate(permutations(1:n))
        if j < noprintterms || j > numperm-noprintterms
            s = Combinatorics.parity(p)
            (j > 1 || s == 1) && write(io, " ", s == 0 ? "+" : "-", " ")
            for (i,c) in enumerate(p)
                write(io, "$(sd.orbitals[i])($(c))")
            end
        elseif j == noprintterms
            write(io, " + … ")
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

export SlaterDeterminant
