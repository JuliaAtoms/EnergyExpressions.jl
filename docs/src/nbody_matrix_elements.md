# N-body matrix elements

```@meta
CurrentModule = EnergyExpressions
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```

The matrix element of an N-body operator between two Slater
determinants may be expanded according to the Löwdin rules (which
reduce to the Slater–Condon rules if all single-particle orbitals are
orthogonal):

$$\begin{equation}
\matrixel{\Phi_A}{\Omega_n}{\Phi_B} =
\sum_p (-)^p
\matrixel{k_1k_2...k_n}{\Omega_n}{l_1l_2...l_n}
D^{AB}({k_1k_2...k_n}|{l_1l_2...l_n})
\end{equation}$$

where $D^{AB}({k_1k_2...k_n}|{l_1l_2...l_n})$ is the determinant minor
of the orbital overlap determinant $D^{AB}$ with the rows
${k_1k_2...k_n}$ and columns ${l_1l_2...l_n}$ stricken out, and $p$
runs over all permutations.

In general, a term in the expansion is thus of the form

$$\begin{equation}
c\matrixel{k_1k_2...k_n}{\Omega_n}{l_1l_2...l_n}\braket{a}{b}\braket{c}{d}\dots\braket{y}{z},
\end{equation}$$

where $c$ is a scalar. This is represented by [`NBodyTerm`](@ref) type.

```@docs
NBodyTermFactor
OrbitalOverlap
OrbitalMatrixElement
NBodyTerm
NBodyMatrixElement
isabovediagonal
isdiagonal
detaxis
detminor
# indexsum
cofactor
det
permutation_sign
overlap_matrix
Matrix
```

```@meta
 DocTestSetup = nothing
```
