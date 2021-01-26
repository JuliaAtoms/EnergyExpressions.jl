# * Common operators

"""
    OneBodyHamiltonian

The one-body Hamiltonian, may include external fields. It is diagonal
in spin, i.e. it does not couple orbitals of opposite spin.
"""
struct OneBodyHamiltonian <: OneBodyOperator end

Base.show(io::IO, ::OneBodyHamiltonian) = write(io, "ĥ")
Base.show(io::IO, me::OrbitalMatrixElement{1,A,OneBodyHamiltonian,B}) where{A,B} =
    write(io, "(", join(string.(me.a), " "), "|", join(string.(me.b), " "), ")")

"""
    iszero(me::EnergyExpressions.OrbitalMatrixElement{1,<:SpinOrbital{<:Orbital},OneBodyHamiltonian,<:SpinOrbital{<:Orbital}})

The matrix element vanishes if the spin-orbitals do not have the same spin.
"""
Base.iszero(me::OrbitalMatrixElement{1,A,OneBodyHamiltonian,B}) where {A<:SpinOrbital{<:Orbital},B<:SpinOrbital{<:Orbital}} =
    me.a[1].m[2] != me.b[1].m[2]

"""
    FieldFreeOneBodyHamiltonian

The one-body Hamiltonian, with no external fields. It is diagonal in
the orbital symmetry.
"""
struct FieldFreeOneBodyHamiltonian <: OneBodyOperator end

Base.show(io::IO, ::FieldFreeOneBodyHamiltonian) = write(io, "ĥ₀")
Base.iszero(me::OrbitalMatrixElement{1,<:SpinOrbital,FieldFreeOneBodyHamiltonian,<:SpinOrbital}) =
    symmetry(me.a[1]) != symmetry(me.b[1])

"""
    CoulombInteraction

Two-body Hamiltonian, representing the mutual Coulombic repulsion
between two electrons. Is diagonal in spin, i.e. the spin of the
orbitals associated with the same coordinate must be the same.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))
[a b|c d]

julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:b,:a))
G(a,b)
```
"""
struct CoulombInteraction{O} <: TwoBodyOperator
    o::O # This can be used to indicate e.g. approximate Coulomb interaction
end
CoulombInteraction() = CoulombInteraction(nothing)

Base.hash(g::CoulombInteraction, h::UInt) = hash(g.o, h)

const CoulombPotential{A,B} = ContractedOperator{1,2,1,A,<:CoulombInteraction,B}

Base.show(io::IO, ::CoulombInteraction) = write(io, "ĝ")

function Base.show(io::IO, me::OrbitalMatrixElement{2,A,<:CoulombInteraction,B}) where {A,B}
    if me.a == me.b # Direct interaction
        write(io, "F($(me.a[1]),$(me.a[2]))")
    elseif me.a[1] == me.b[2] && me.a[2] == me.b[1] # Exchange interaction
        write(io, "G($(me.a[1]),$(me.a[2]))")
    else # General case
        write(io, "[", join(string.(me.a), " "), "|", join(string.(me.b), " "), "]")
    end
end

Base.show(io::IO, co::CoulombPotential{A,B}) where {A,B} =
    write(io, "[$(co.a[1])|$(co.b[1])]")

"""
    iszero(me::EnergyExpressions.OrbitalMatrixElement{2,<:SpinOrbital{<:Orbital},CoulombInteraction,<:SpinOrbital{<:Orbital}})

The matrix element vanishes if the (non-relativistic) spin-orbitals
associated with the same coordinate do not have the same spin.
"""
Base.iszero(me::OrbitalMatrixElement{2,A,<:CoulombInteraction,B}) where {A<:SpinOrbital{<:Orbital},B<:SpinOrbital{<:Orbital}} =
    me.a[1].m[2] != me.b[1].m[2] || me.a[2].m[2] != me.b[2].m[2]

export OneBodyHamiltonian, FieldFreeOneBodyHamiltonian, CoulombInteraction, CoulombPotential
