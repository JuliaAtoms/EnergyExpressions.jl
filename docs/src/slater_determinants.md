# Slater determinants

[Slater
determinants](https://en.wikipedia.org/wiki/Slater_determinant) are
wavefunctions constructed from anti-symmetrized one-particle states.

```@meta
CurrentModule = EnergyExpressions
```

## Construction of Slater determinants
```@docs
SlaterDeterminant
SlaterDeterminant(::Configuration{<:SpinOrbital})
length(::SlaterDeterminant)
AdjointSlaterDeterminant
adjoint(::SlaterDeterminant)
length(::AdjointSlaterDeterminant)
```

## N-body operators
```@docs
NBodyOperator
NBodyTermFactor
OrbitalOverlap
OrbitalMatrixElement
```

## Example usage

```jldoctest
julia> sa = SlaterDeterminant([:l, :a])
l(1)a(2) - l(2)a(1)

julia> sb = SlaterDeterminant([:k, :b])
k(1)b(2) - k(2)b(1)

julia> sa'ZeroBodyOperator*sb
+ ⟨l|k⟩⟨a|b⟩ - ⟨l|b⟩⟨a|k⟩ - ⟨a|k⟩⟨l|b⟩ + ⟨a|b⟩⟨l|k⟩

julia> sa'OneBodyOperator*sb
+ ⟨l|Ω₁|k⟩⟨a|b⟩ + ⟨l|k⟩⟨a|Ω₁|b⟩ - ⟨l|Ω₁|b⟩⟨a|k⟩ - ⟨l|b⟩⟨a|Ω₁|k⟩ - ⟨a|Ω₁|k⟩⟨l|b⟩ - ⟨a|k⟩⟨l|Ω₁|b⟩ + ⟨a|Ω₁|b⟩⟨l|k⟩ + ⟨a|b⟩⟨l|Ω₁|k⟩
```
