# N-body equations

```@meta
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```
These types are used to represent the integro-differential equations
that result from the variation of the energy expression with respect
to an orbital. They can, in general, be written on the form

$$\begin{equation}
0 = \Omega_1\ket{\chi}\braket{a}{b}...
\end{equation}$$

when varying with respect to a conjugated orbital, and

$$\begin{equation}
0 = \bra{\chi}\Omega_1\braket{a}{b}...
\end{equation}$$

when varying with respect to an orbital. In both cases, $\Omega_1$ is
a one-body operator, either in itself, or resulting from a contraction
over all coordinates but one of a many-body operator (see
[`ContractedOperator`](@ref)).

```@docs
NBodyEquation
LinearCombinationEquation
```

```@meta
 DocTestSetup = nothing
```
