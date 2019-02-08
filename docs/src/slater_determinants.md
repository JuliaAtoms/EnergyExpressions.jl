# Slater determinants

[Slater
determinants](https://en.wikipedia.org/wiki/Slater_determinant) are
wavefunctions constructed from anti-symmetrized one-particle states.

```@meta
CurrentModule = EnergyExpressions
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
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

## Example usage

```jldoctest
julia> sa = SlaterDeterminant([:l, :a])
l(1)a(2) - l(2)a(1)

julia> sb = SlaterDeterminant([:k, :b])
k(1)b(2) - k(2)b(1)

julia> using AtomicLevels

julia> SlaterDeterminant(spin_configurations(c"1s2")[1])
1s₀α(1)1s₀β(2) - 1s₀α(2)1s₀β(1)
```

```@meta
 DocTestSetup = nothing
```
