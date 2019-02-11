# System of equations

When dealing with large energy expressions resulting from many
configurations, when deriving the orbital equations
(see [`N-body equations`](@ref) and
[`Calculus of Variations – Implementation`](@ref)), many of the
integrals involved will be shared between the different terms of the
equations. It is therefore of interest to gather a list of integrals
that can be calculated once at every iteration (in a self-consistent
procedure for finding eigenstates, or in time-propagation), and whose
values can then be reused when solving the individual equations. The
routines described below, although primitive, aid in this effort.

The idea is the following: With an energy expression on the form

$$\begin{equation}
E(\vec{P},\vec{c}) = \frac{\vec{c}^H
\mat{H}\vec{c}}{\vec{c}^H\vec{c}}
\end{equation}$$

where $\vec{P}$ is a set of orbitals, $\vec{c}$ is a vector of
mixing coefficients, and

$$\begin{equation}
\mat{H}_{ij} \defd \matrixel{\chi_i}{\Hamiltonian}{\chi_j},
\end{equation}$$

all terms of the equations can be written on the form

$$\begin{equation}
\label{eqn:equation-term}
\operator{A}\ket{\chi}
\underbrace{\left[
\sum_n \alpha_n \conj{c}_{i_n} c_{j_n}
\prod_k \int_{n_k}
\right]}_{\textrm{Rank 0}},
\end{equation}$$

(or on its dual form with conjugated orbitals) where $\operator{A}$ is
some one-body operator, $\alpha_n$ a coefficient due to the energy
expression, $\conj{c}_{i_n} c_{j_n}$ the coefficient due to the
multi-configurational expansion, and finally $\int_{n_k}$ all the
extra integrals (which are necessarily zero-body operators) that may
appear due to non-orthogonal orbitals (and which may be shared between
many equations). The operator $\operator{A}$ can have ranks 0–2;
formally, it can only have rank 0 (multiplication by a scalar) or 2
(multiplication by a matrix). What we somewhat sloppily to refer by
“rank 1” is an operator of rank 2, but which is diagonal in the
underlying coordinate, such as a local potential.

Thus, to efficiently perform an iteration, one would first compute all
common integrals $\int_k$, and then for every equation term of the
form $\eqref{eqn:equation-term}$, form the coefficient in the brakets,
calculate the action of the operator $\operator{A}$ on the orbital
$\chi$, multiplied the coefficient.

There are a few important special cases of
$\eqref{eqn:equation-term}$:

1) Field-free, one-body Hamiltonian, i.e. $\operator{A}=\hamiltonian_0$.
   This is the contribution from the orbital $\ket{\chi}$ to
   itself. Rank 2.

2) One-body Hamiltonian, including field, i.e. $\operator{A}=\hamiltonian$.
   This is the contribution from another orbital $\ket{\chi'}$ to
   $\ket{\chi}$ via some off-diagonal coupling, such as an external
   field (e.g. an electro–magnetic field). Rank 2.

3) Direct interaction, i.e. $\operator{A}=\direct{kl}$, where two
   other orbitals $\ket{\chi_k}$ and $\ket{\chi_l}$ together form a
   potential acting on $\ket{\chi}$. “Rank 1”.

4) Exchange interaction, i.e. $\operator{A}=\exchange{kl}$, where
   another orbital $\ket{\chi_k}$ and $\ket{\chi}$ together form a
   potential acting on a third orbital $\ket{\chi_l}$. Rank 2.

5) Source term, i.e. a contribution that does not involve $\ket{\chi}$
   in any way. This term arises from other configurations in the
   multi-configurational expansion. Case 2. is also formulated in this
   way. Rank 0–2. If for some reason the source orbital and/or the
   operator acting on it is fixed, it may be possible to precompute
   the effect of the operator on the source orbital and reduce the
   computational complexity to a rank 0-equivalent operation.

```@meta
CurrentModule = EnergyExpressions
DocTestSetup = quote
    using EnergyExpressions
end
```

## Example

We start by defining a set of configurations and specifying that a few
of the constituent orbitals are non-orthogonal:

```jldoctest mc-eqs
julia> cfgs = [[:i, :j, :l, :k̃], [:i, :j, :k, :l̃]]
2-element Array{Array{Symbol,1},1}:
 [:i, :j, :l, :k̃]
 [:i, :j, :k, :l̃]

julia> continua = [:k̃, :l̃]
2-element Array{Symbol,1}:
 :k̃
 :l̃

julia> overlaps = [OrbitalOverlap(i,j) for i in continua for j in continua]
4-element Array{OrbitalOverlap{Symbol,Symbol},1}:
 ⟨k̃|k̃⟩
 ⟨k̃|l̃⟩
 ⟨l̃|k̃⟩
 ⟨l̃|l̃⟩
```

We then set up the energy expression as before:

```jldoctest mc-eqs
julia> H = OneBodyHamiltonian() + CoulombInteraction()
ĥ + ĝ

julia> E = Matrix(H, SlaterDeterminant.(cfgs), overlaps)
2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (i|i)⟨k̃|k̃⟩ + (j|j)⟨k̃|k̃⟩ + (l|l)⟨k̃|k̃⟩ + (k̃|k̃) - G(i,j)⟨k̃|k̃⟩ + F(i,j)⟨k̃|k̃⟩ + … - G(i,k̃) + F(i,k̃) - G(j,k̃) + F(j,k̃) - G(l,k̃) + F(l,k̃)  …  (l|k)⟨k̃|l̃⟩ - [i l|k i]⟨k̃|l̃⟩ + [i l|i k]⟨k̃|l̃⟩ - [j l|k j]⟨k̃|l̃⟩ + [j l|j k]⟨k̃|l̃⟩ - [l k̃|l̃ k] + [l k̃|k l̃]
 (k|l)⟨l̃|k̃⟩ - [i k|l i]⟨l̃|k̃⟩ + [i k|i l]⟨l̃|k̃⟩ - [j k|l j]⟨l̃|k̃⟩ + [j k|j l]⟨l̃|k̃⟩ - [k l̃|k̃ l] + [k l̃|l k̃]                                     (i|i)⟨l̃|l̃⟩ + (j|j)⟨l̃|l̃⟩ + (k|k)⟨l̃|l̃⟩ + (l̃|l̃) - G(i,j)⟨l̃|l̃⟩ + F(i,j)⟨l̃|l̃⟩ + … - G(i,l̃) + F(i,l̃) - G(j,l̃) + F(j,l̃) - G(k,l̃) + F(k,l̃)
```

Finally, we derive the coupled integro-differential equation system
for the continuum orbitals `k̃`, `l̃`:

```jldoctest mc-eqs
julia> eqs = diff(E, Conjugate.(continua))
EnergyExpressions.MCEquationSystem(EnergyExpressions.OrbitalEquation{Symbol,SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64}}[OrbitalEquation(k̃):
  [1, 1]  =  𝐈₁|k̃⟩(i|i) + 𝐈₁|k̃⟩(j|j) + 𝐈₁|k̃⟩(l|l) + ĥ|k̃⟩ + 𝐈₁|k̃⟩(- G(i,j)) + 𝐈₁|k̃⟩F(i,j) + 𝐈₁|k̃⟩(- G(i,l)) + 𝐈₁|k̃⟩F(i,l) + 𝐈₁|k̃⟩(- G(j,l)) + 𝐈₁|k̃⟩F(j,l) + -[i|k̃]|i⟩ + [i|i]|k̃⟩ + -[j|k̃]|j⟩ + [j|j]|k̃⟩ + -[l|k̃]|l⟩ + [l|l]|k̃⟩
  [1, 2]  =  𝐈₁|l̃⟩(l|k) + 𝐈₁|l̃⟩(- [i l|k i]) + 𝐈₁|l̃⟩[i l|i k] + 𝐈₁|l̃⟩(- [j l|k j]) + 𝐈₁|l̃⟩[j l|j k] + -[l|l̃]|k⟩ + [l|k]|l̃⟩
, OrbitalEquation(l̃):
  [2, 1]  =  𝐈₁|k̃⟩(k|l) + 𝐈₁|k̃⟩(- [i k|l i]) + 𝐈₁|k̃⟩[i k|i l] + 𝐈₁|k̃⟩(- [j k|l j]) + 𝐈₁|k̃⟩[j k|j l] + -[k|k̃]|l⟩ + [k|l]|k̃⟩
  [2, 2]  =  𝐈₁|l̃⟩(i|i) + 𝐈₁|l̃⟩(j|j) + 𝐈₁|l̃⟩(k|k) + ĥ|l̃⟩ + 𝐈₁|l̃⟩(- G(i,j)) + 𝐈₁|l̃⟩F(i,j) + 𝐈₁|l̃⟩(- G(i,k)) + 𝐈₁|l̃⟩F(i,k) + 𝐈₁|l̃⟩(- G(j,k)) + 𝐈₁|l̃⟩F(j,k) + -[i|l̃]|i⟩ + [i|i]|l̃⟩ + -[j|l̃]|j⟩ + [j|j]|l̃⟩ + -[k|l̃]|k⟩ + [k|k]|l̃⟩
], Any[(i|i), 𝐈₁, (j|j), (l|l), ĥ, G(i,j), F(i,j), G(i,l), F(i,l), G(j,l)  …  [j k|l j], [j k|j l], [k|k̃], [k|l], (k|k), G(i,k), F(i,k), G(j,k), F(j,k), [k|k]])
```

We can investigate the [`MCEquationSystem`](@ref) object `eqs` a
bit. It consists of two coupled equations:

```jldoctest mc-eqs
julia> eqs.equations
2-element Array{EnergyExpressions.OrbitalEquation{Symbol,SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64}},1}:
 OrbitalEquation(k̃):
  [1, 1]  =  𝐈₁|k̃⟩(i|i) + 𝐈₁|k̃⟩(j|j) + 𝐈₁|k̃⟩(l|l) + ĥ|k̃⟩ + 𝐈₁|k̃⟩(- G(i,j)) + 𝐈₁|k̃⟩F(i,j) + 𝐈₁|k̃⟩(- G(i,l)) + 𝐈₁|k̃⟩F(i,l) + 𝐈₁|k̃⟩(- G(j,l)) + 𝐈₁|k̃⟩F(j,l) + -[i|k̃]|i⟩ + [i|i]|k̃⟩ + -[j|k̃]|j⟩ + [j|j]|k̃⟩ + -[l|k̃]|l⟩ + [l|l]|k̃⟩
  [1, 2]  =  𝐈₁|l̃⟩(l|k) + 𝐈₁|l̃⟩(- [i l|k i]) + 𝐈₁|l̃⟩[i l|i k] + 𝐈₁|l̃⟩(- [j l|k j]) + 𝐈₁|l̃⟩[j l|j k] + -[l|l̃]|k⟩ + [l|k]|l̃⟩

 OrbitalEquation(l̃):
  [2, 1]  =  𝐈₁|k̃⟩(k|l) + 𝐈₁|k̃⟩(- [i k|l i]) + 𝐈₁|k̃⟩[i k|i l] + 𝐈₁|k̃⟩(- [j k|l j]) + 𝐈₁|k̃⟩[j k|j l] + -[k|k̃]|l⟩ + [k|l]|k̃⟩
  [2, 2]  =  𝐈₁|l̃⟩(i|i) + 𝐈₁|l̃⟩(j|j) + 𝐈₁|l̃⟩(k|k) + ĥ|l̃⟩ + 𝐈₁|l̃⟩(- G(i,j)) + 𝐈₁|l̃⟩F(i,j) + 𝐈₁|l̃⟩(- G(i,k)) + 𝐈₁|l̃⟩F(i,k) + 𝐈₁|l̃⟩(- G(j,k)) + 𝐈₁|l̃⟩F(j,k) + -[i|l̃]|i⟩ + [i|i]|l̃⟩ + -[j|l̃]|j⟩ + [j|j]|l̃⟩ + -[k|l̃]|k⟩ + [k|k]|l̃⟩
```

The first equation consists of the following terms:

```jldoctest mc-eqs
julia> eqs.equations[1].one_body
0-element Array{EnergyExpressions.MCCoeff,1}

julia> eqs.equations[1].direct_terms
3-element Array{Pair{Int64,Array{EnergyExpressions.MCCoeff,1}},1}:
 13 => [MCCoeff{Int64}(1, 1, 1, Int64[])]
 14 => [MCCoeff{Int64}(1, 1, 1, Int64[])]
 12 => [MCCoeff{Int64}(1, 1, 1, Int64[])]

julia> eqs.equations[1].exchange_terms
3-element Array{Tuple{Any,EnergyExpressions.MCCoeff,Any},1}:
 ([i|k̃], EnergyExpressions.MCCoeff{Int64}(1, 1, -1, Int64[]), :i)
 ([j|k̃], EnergyExpressions.MCCoeff{Int64}(1, 1, -1, Int64[]), :j)
 ([l|k̃], EnergyExpressions.MCCoeff{Int64}(1, 1, -1, Int64[]), :l)

julia> eqs.equations[1].source_terms
4-element Array{Pair{Int64,Array{Tuple{EnergyExpressions.MCCoeff,Any},1}},1}:
  2 => [(MCCoeff{Int64}(1, 1, 1, [1]), :k̃), (MCCoeff{Int64}(1, 1, 1, [3]), :k̃), (MCCoeff{Int64}(1, 1, 1, [4]), :k̃), (MCCoeff{Int64}(1, 1, -1, [6]), :k̃), (MCCoeff{Int64}(1, 1, 1, [7]), :k̃), (MCCoeff{Int64}(1, 1, -1, [8]), :k̃), (MCCoeff{Int64}(1, 1, 1, [9]), :k̃), (MCCoeff{Int64}(1, 1, -1, [10]), :k̃), (MCCoeff{Int64}(1, 1, 1, [11]), :k̃), (MCCoeff{Int64}(1, 2, 1, [15]), :l̃), (MCCoeff{Int64}(1, 2, -1, [16]), :l̃), (MCCoeff{Int64}(1, 2, 1, [17]), :l̃), (MCCoeff{Int64}(1, 2, -1, [18]), :l̃), (MCCoeff{Int64}(1, 2, 1, [19]), :l̃)]
  5 => [(MCCoeff{Int64}(1, 1, 1, Int64[]), :k̃)]
 21 => [(MCCoeff{Int64}(1, 2, 1, Int64[]), :l̃)]
 20 => [(MCCoeff{Int64}(1, 2, -1, Int64[]), :k)]
```

where the [`MCCoeff`](@ref) objects indicate which components of the
mixing coefficient vector $\vec{c}$ need to be multiplied, and all the
integers are pointers to the list of common integrals:

```jldoctest mc-eqs
julia> eqs.integrals
34-element Array{Any,1}:
 (i|i)
 𝐈₁
 (j|j)
 (l|l)
 ĥ
 G(i,j)
 F(i,j)
 G(i,l)
 F(i,l)
 G(j,l)
 F(j,l)
 [i|i]
 [j|j]
 [l|l]
 (l|k)
 [i l|k i]
 [i l|i k]
 [j l|k j]
 [j l|j k]
 [l|l̃]
 [l|k]
 (k|l)
 [i k|l i]
 [i k|i l]
 [j k|l j]
 [j k|j l]
 [k|k̃]
 [k|l]
 (k|k)
 G(i,k)
 F(i,k)
 G(j,k)
 F(j,k)
 [k|k]
```

From this we see that the `𝐈₁` (one-body identity operator)
contribution to `|k̃⟩` can be written as a linear combination of `|k̃⟩`
and `|l̃⟩`, weighted by different components of the mixing coefficient
vector $\vec{c}$ and various other integrals. This is all the
information necessary to set up an efficient equation solver.

## Implementation

```@docs
MCCoeff
OrbitalEquation
MCEquationSystem
pushifmissing!
```

```@meta
 DocTestSetup = nothing
```

