# N-body matrix elements

```@meta
CurrentModule = EnergyExpressions
DocTestSetup = quote
    using EnergyExpressions
    using AtomicLevels
end
```

The matrix element of an $N$-body operator between two Slater
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
numbodies
NBodyTerm
NBodyMatrixElement
isdependent
transform
overlap_matrix
EnergyExpressions.compare
Base.:(==)(::EnergyExpressions.NBodyMatrixElement, ::EnergyExpressions.NBodyMatrixElement)
Base.isapprox(::EnergyExpressions.NBodyMatrixElement, ::EnergyExpressions.NBodyMatrixElement)
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
permutations, as well as the fact that the overlap matrix is very
sparse.

### Finding non-zero minors of the overlap determinant

The algorithm to find which minor determinants
$\Gamma^{(N)}(k_1k_2...k_N|l_1l_2...l_N)$ do not vanish, and hence
which $N$ orbitals $k_1k_2...k_N,l_1l_2...l_N$ the $N$-body operator
should be contracted over, is described briefly below. It is devised
to be optimal for orthogonal orbitals (i.e. linear complexity
$\mathcal{O}(Nn)$ where $n$ is the number of orbitals), and
near-optimal for a small amount of non-orthogonal orbitals.

Given:
- An $N$-body operator (implying $N$ rows and $N$ need to
  stricken out), and
- an $n\times n$ matrix, with coordinates of non-zero matrix
  elements: $I,J$ (from these vectors, vanishing rows/columns can
  easily be deduced $\implies$ $N_r$ row/$N_c$ column "rank",
  i.e. yet to be stricken out),
do
- find all $N_r$-combinations of the remaining rows,
- for each such combination, find all columns which would be
  affected if striking out that particular combination of rows,
  - if the "support" (i.e. the only non-zero elements) of any of those
    columns vanishes when striking out the rows, that column must be
    stricken out, too. Total number of these columns is named
    $N_{cm}$,
  - if more than $N_c$ columns must be stricken out ($N_{cm}>N_c$), that
    row combination is unviable,
  - find all $N_c - N_{cm}$-combinations of the remaining columns,
  - for each such combination of columns, find all rows which would be
    affected,
    - if the support of any of those rows vanishes when striking out the
      candidate columns, and the row is not in the candidate set of rows
      to be stricken out, the column combination is unviable.

```@docs
nonzero_minors
```

### Contracting over the stricken out orbitals

We use Julia's built-in `Base.Cartesian.@nloops` iterators to span the
space of all possible choices of orbitals for the contraction of
orbitals. If two or more orbitals are the same, the matrix element is
trivially zero (the “Fermi hole”). To avoid double-counting, we also
only consider those indices that are above the hyper-diagonal.

```@docs
detaxis
detminor
# indexsum
cofactor
EnergyExpressions.distinct_permutations
det
permutation_sign
powneg1
EnergyExpressions.@above_diagonal_loop
EnergyExpressions.@anti_diagonal_loop
```


## Occupation number representation

As an alternative to [`SlaterDeterminant`](@ref), we also provide an
implementation of the occupation number representation (i.e. second
quantization). Every configuration is thus represented as a bit
vector, indicating the occupation of a specific orbital by a `true`
value. Excitations between two configurations are easily computed
using bitwise operations. The present implementation is inspired by

- Scemama, A., & Giner, E. (2013). An efficient implementation of
  Slater–Condon rules. CoRR,
  [arXiv:1311.6244](https://arxiv.org/abs/1311.6244),

but extends upon it by also supporting non-orthogonal orbitals.

The benefit of this approach is much more efficient identification of
which cofactors in $\eqref{eqn:matrix-element-expansion}$ are
non-zero, than the approach taken in [`nonzero_minors`](@ref).

```@docs
BitConfigurations
EnergyExpressions.Orbitals
EnergyExpressions.non_zero_cofactors
```

```@meta
 DocTestSetup = nothing
```
