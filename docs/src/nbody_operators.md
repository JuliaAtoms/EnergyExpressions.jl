# N-body operators

```@meta
CurrentModule = EnergyExpressions
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```

```@docs
NBodyOperator
LinearCombinationOperator
IdentityOperator
ContractedOperator
Base.in(orbital::O, co::ContractedOperator) where O
Base.in(corbital::Conjugate{O}, co::ContractedOperator) where O
contract
complement
LinearAlgebra.ishermitian(::QuantumOperator)
```

```@meta
 DocTestSetup = nothing
```
