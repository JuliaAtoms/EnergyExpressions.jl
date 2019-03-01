# Conjugate orbitals

```@meta
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```

A conjugated orbital, $\conj{\chi}$ (often written $\bra{\chi}$) is
the dual to the unconjugated orbital $\chi$ (often written
$\ket{\chi}$). In the code, conjugation of orbitals is denoted with a
dagger (`†`) to avoid confusion with the multiplication operator `*`.

```@docs
Conjugate
conj
```

```@meta
 DocTestSetup = nothing
```
