# * Common operators

"""
    OneBodyHamiltonian

The one-body Hamiltonian, may include external fields. Is diagonal in spin, i.e.
"""
struct OneBodyHamiltonian <: OneBodyOperator end

Base.show(io::IO, ::OneBodyHamiltonian) = write(io, "ĥ")
Base.show(io::IO, me::EnergyExpressions.OrbitalMatrixElement{1,A,OneBodyHamiltonian,B}) where{A,B} =
    write(io, "(", join(string.(me.a), " "), "|", join(string.(me.b), " "), ")")

"""
    iszero(me::EnergyExpressions.OrbitalMatrixElement{1,<:SpinOrbital,OneBodyHamiltonian,<:SpinOrbital})

The matrix element vanishes if the spin-orbitals do not have the same spin.
"""
Base.iszero(me::EnergyExpressions.OrbitalMatrixElement{1,A,OneBodyHamiltonian,B}) where {A<:SpinOrbital,B<:SpinOrbital} =
    me.a[1].spin != me.b[1].spin

"""
    FieldFreeOneBodyHamiltonian

The one-body Hamiltonian, with no external fields. Is diagonal in the
orbitals, i.e. does not couple unequal orbitals.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))
[a b|c d]
```
"""
struct FieldFreeOneBodyHamiltonian <: OneBodyOperator end

Base.show(io::IO, ::FieldFreeOneBodyHamiltonian) = write(io, "ĥ₀")
Base.iszero(me::EnergyExpressions.OrbitalMatrixElement{1,A,FieldFreeOneBodyHamiltonian,B}) where {A,B} =
    me.a != me.b

"""
    CoulombInteraction

Two-body Hamiltonian, representing the mutual Coulombic repulsion
between two electrons. Is diagonal in spin, i.e. the spin of the
orbitals associated with the same coordinate must be the same.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))
[a b|c d]
```
"""
struct CoulombInteraction <: TwoBodyOperator end

Base.show(io::IO, ::CoulombInteraction) = write(io, "ĝ")
Base.show(io::IO, me::EnergyExpressions.OrbitalMatrixElement{2,A,CoulombInteraction,B}) where{A,B} =
    write(io, "[", join(string.(me.a), " "), "|", join(string.(me.b), " "), "]")

"""
    iszero(me::EnergyExpressions.OrbitalMatrixElement{2,<:SpinOrbital,CoulombInteraction,<:SpinOrbital})

The matrix element vanishes if the spin-orbitals associated with the
same coordinate do not have the same spin.
"""
Base.iszero(me::EnergyExpressions.OrbitalMatrixElement{2,A,CoulombInteraction,B}) where {A<:SpinOrbital,B<:SpinOrbital} =
    me.a[1].spin != me.b[1].spin || me.a[2].spin != me.b[2].spin

export OneBodyHamiltonian, FieldFreeOneBodyHamiltonian, CoulombInteraction
