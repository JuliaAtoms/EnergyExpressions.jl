# Common N-body operators

```@meta
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```

A few N-body operators that are common in quantum mechanics. It is
possible to form a composite Hamiltonian thus:

```jldoctest common_operators
julia> H = OneBodyHamiltonian() + CoulombInteraction()
ĥ + ĝ
```
If then define a set of Slater determinants, we can easily form the energy expression:
```jldoctest common_operators
julia> slaters = SlaterDeterminant.([[:a, :b], [:c, :d], [:a, :c]])
3-element Array{SlaterDeterminant{Symbol},1}:
 |a b|
 |c d|
 |a c|

julia> Matrix(H, slaters)
3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (a|a) + (b|b) - G(a,b) + F(a,b)  - [a b|d c] + [a b|c d]          (b|c) - [a b|c a] + [a b|a c]
 - [c d|b a] + [c d|a b]          (c|c) + (d|d) - G(c,d) + F(c,d)  - (d|a) - [c d|c a] + [c d|a c]
 (c|b) - [a c|b a] + [a c|a b]    - (a|d) - [a c|d c] + [a c|c d]  (a|a) + (c|c) - G(a,c) + F(a,c)
```
An energy expression like this can then be used to derive the
multi-configurational Hartree–Fock equations for the orbitals
`a,b,c,d`.

We can also specify that e.g. the orbitals `b` and `c` are
non-orthogonal, and thus derive a slightly different energy
expression, that takes this into account:
```jldoctest common_operators
julia> Matrix(H, slaters, [OrbitalOverlap(:b,:c)])
3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (a|a) + (b|b) - G(a,b) + F(a,b)             - ⟨b|c⟩(a|d) - [a b|d c] + [a b|c d]  ⟨b|c⟩(a|a) + (b|c) - [a b|c a] + [a b|a c]
 - ⟨c|b⟩(d|a) - [c d|b a] + [c d|a b]        (c|c) + (d|d) - G(c,d) + F(c,d)       - (d|a) - [c d|c a] + [c d|a c]
 ⟨c|b⟩(a|a) + (c|b) - [a c|b a] + [a c|a b]  - (a|d) - [a c|d c] + [a c|c d]       (a|a) + (c|c) - G(a,c) + F(a,c)
```
Beware that the computational complexity grows factorially with the
amount of non-orthogonal orbitals!

```@docs
OneBodyHamiltonian
FieldFreeOneBodyHamiltonian
CoulombInteraction
iszero
```

```@meta
 DocTestSetup = nothing
```
