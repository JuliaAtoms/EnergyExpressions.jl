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
\label{eqn:matrix-element-expansion}
\matrixel{\Phi_A}{\Omega_n}{\Phi_B} =
\frac{1}{n!}\sum_p (-)^p
\matrixel{k_1k_2...k_n}{\Omega_n}{l_1l_2...l_n}
D^{AB}({k_1k_2...k_n}|{l_1l_2...l_n})
\end{equation}$$

where $D^{AB}({k_1k_2...k_n}|{l_1l_2...l_n})$ is the determinant minor
of the orbital overlap determinant $D^{AB}$ with the rows
${k_1k_2...k_n}$ and columns ${l_1l_2...l_n}$ stricken out, and $p$
runs over all permutations.

In general, a term in the expansion is thus of the form

$$\begin{equation}
\alpha\matrixel{k_1k_2...k_n}{\Omega_n}{l_1l_2...l_n}\braket{a}{b}\braket{c}{d}\dots\braket{y}{z},
\end{equation}$$

where $\alpha$ is a scalar. This is represented by [`NBodyTerm`](@ref) type.

```@docs
NBodyTermFactor
OrbitalOverlap
OrbitalMatrixElement
NBodyTerm
NBodyMatrixElement
isdependent
transform
overlap_matrix
EnergyExpression
Matrix
```

## Calculation of determinants

Actually computing the matrix element expansion
$\eqref{eqn:matrix-element-expansion}$ is a combinatorial problem,
that grows factorially with the amount of non-orthogonal orbital
pairs. Furthermore, of the $(n!)^2$ terms generated from the
expansion, only $n!$ are distinct, due to the integrals being
symmetric with respect to interchange of the coordinates [hence the
normalization factor $(n!)^{-1}$]. Thankfully, there are few
symmetries that can be employed, to generate only the distinct
permutations.

We use Julia's built in `CartesianIndex` iterators to span the space
of all possible choices of orbitals for the overlap determinant. If
two or more indices in the `CartesianIndex` are the same, the overlap
is trivially zero (the “Fermi hole”). To avoid double-counting, we
also only consider those indices that are above the hyper-diagonal.

```@docs
isabovediagonal
isdiagonal
detaxis
detminor
# indexsum
cofactor
det
permutation_sign
```

```@meta
 DocTestSetup = nothing
```
