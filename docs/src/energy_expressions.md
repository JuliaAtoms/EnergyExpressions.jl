# Energy expressions

The average energy of the system is given by

$$\begin{equation}
E_{\textrm{av}} \defd \matrixel{\Psi}{\Hamiltonian}{\Psi},
\end{equation}$$

where $\Psi$ is the (multi-electron) wavefunction of the system and
$\Hamiltonian$ is the full Hamiltonian. A common approach is to
approximate the wavefunction as a linear combination of [Slater
determinants](https://en.wikipedia.org/wiki/Slater_determinant):

$$\begin{equation}
\Psi \approx \sum_K D_K\Phi_K,
\end{equation}$$

where each $\Phi_K$ constitutes an anti-symmetrized product of
one-particle spin-orbitals $\chi_i$. Depending on whether these
spin-orbitals are chosen to be orthogonal or not, with respect to each
other, the energy expression takes different forms.

Other choices for the expansion are possible, for instance
[Configuration state
functions](https://en.wikipedia.org/wiki/Configuration_state_function). The
algebra for generating energy expressions from such objects is more
involved, and is not implemented in this library.

The examples below are for [atomic
configurations](https://github.com/JuliaAtoms/AtomicLevels.jl), but
the library is not limited to this case, rather, it could work with
any kind of configuration constructed from single-particle
spin-orbitals. For example, it possible to work with normal Julia
`Symbol`s denoting the spin-orbitals. This can be used to derive
equations of motion purely symbolically.

## Orthogonal case, Slater–Condon rules

```@meta
DocTestSetup = quote
    using AtomicLevels
    using EnergyExpressions
end
```

First, we find all configurations possible for two configurations of
helium, and convert them into [`SlaterDeterminant`](@ref)s.

```jldoctest helium
julia> he = SlaterDeterminant.(spin_configurations(c"1s2"))
1-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:
 |1s₀α 1s₀β|

julia> he_exc = SlaterDeterminant.(spin_configurations(c"1s 2p"))
12-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:
 |1s₀α 2p₋₁α|
 |1s₀α 2p₋₁β|
 |1s₀α 2p₀α|
 |1s₀α 2p₀β|
 |1s₀α 2p₁α|
 |1s₀α 2p₁β|
 |1s₀β 2p₋₁α|
 |1s₀β 2p₋₁β|
 |1s₀β 2p₀α|
 |1s₀β 2p₀β|
 |1s₀β 2p₁α|
 |1s₀β 2p₁β|
```

### One-body energy term

First we create a [`OneBodyHamiltonian`](@ref), the one-body operator
that corresponds to the physical quantity _energy_:

```jldoctest helium
julia> h = OneBodyHamiltonian()
ĥ
```

We can then evaluate matrix elements between Slater determinants
corresponding to different configurations.

The Slater–Condon rules state that the one-body energy expression between two
configurations is

- in case the configurations are identical, $\onebody{\Phi_A}{\Phi_B} = \onebody{i}{i}$:
  ```jldoctest helium
  julia> Matrix(h, he)
  1×1 Array{EnergyExpressions.NBodyMatrixElement,2}:
   (1s₀α|1s₀α) + (1s₀β|1s₀β)
  ```
- in case the configurations differ by one orbital, $\onebody{\Phi_A}{\Phi_B} = \onebody{i}{j}$:
  ```jldoctest helium
  julia> Matrix(h, [he[1],he_exc[2]])
  2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:
   (1s₀α|1s₀α) + (1s₀β|1s₀β)  (1s₀β|2p₋₁β)
   (2p₋₁β|1s₀β)               (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β)
  ```
  We had to choose the second excited state, since the
  [`OneBodyHamiltonian`](@ref) does not couple orbitals of differing
  spin; therefore the corresponding matrix elements would be zero. The
  diagonal matrix elements correspond to the first case.
- in case the configurations differ by more than one orbital, $\onebody{\Phi_A}{\Phi_B} = 0$:
  ```jldoctest helium
  julia> Matrix(h, [he[1],he_exc[10]])
  2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:
   (1s₀α|1s₀α) + (1s₀β|1s₀β)  0
   0                          (1s₀β|1s₀β) + (2p₀β|2p₀β)
  ```

We can easily generate the full one-body matrix:

```jldoctest helium
julia> Matrix(h, vcat(he,he_exc))
13×13 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (1s₀α|1s₀α) + (1s₀β|1s₀β)  0                            (1s₀β|2p₋₁β)                 …  - (1s₀α|2p₁α)              0
 0                          (1s₀α|1s₀α) + (2p₋₁α|2p₋₁α)  0                               0                          0
 (2p₋₁β|1s₀β)               0                            (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β)     0                          0
 0                          (2p₀α|2p₋₁α)                 0                               0                          0
 (2p₀β|1s₀β)                0                            (2p₀β|2p₋₁β)                    0                          0
 0                          (2p₁α|2p₋₁α)                 0                            …  0                          0
 (2p₁β|1s₀β)                0                            (2p₁β|2p₋₁β)                    0                          0
 - (2p₋₁α|1s₀α)             0                            0                               (2p₋₁α|2p₁α)               0
 0                          0                            0                               0                          (2p₋₁β|2p₁β)
 - (2p₀α|1s₀α)              0                            0                               (2p₀α|2p₁α)                0
 0                          0                            0                            …  0                          (2p₀β|2p₁β)
 - (2p₁α|1s₀α)              0                            0                               (1s₀β|1s₀β) + (2p₁α|2p₁α)  0
 0                          0                            0                               0                          (1s₀β|1s₀β) + (2p₁β|2p₁β)
```

### Two-body energy term

Similar considerations apply for the two-body energy terms between two
configurations. To make it more interesting, we consider lithium
which has three electrons:

```jldoctest lithium
julia> li = SlaterDeterminant.(spin_configurations(c"1s2 2s"))
2-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:
 |1s₀α 1s₀β 2s₀α|
 |1s₀α 1s₀β 2s₀β|

julia> li_exc = SlaterDeterminant.(spin_configurations(c"1s 2s 2p"))
24-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:
 |1s₀α 2s₀α 2p₋₁α|
 |1s₀α 2s₀α 2p₋₁β|
 |1s₀α 2s₀α 2p₀α|
 |1s₀α 2s₀α 2p₀β|
 |1s₀α 2s₀α 2p₁α|
 |1s₀α 2s₀α 2p₁β|
 |1s₀α 2s₀β 2p₋₁α|
 |1s₀α 2s₀β 2p₋₁β|
 |1s₀α 2s₀β 2p₀α|
 |1s₀α 2s₀β 2p₀β|
 |1s₀α 2s₀β 2p₁α|
 |1s₀α 2s₀β 2p₁β|
 |1s₀β 2s₀α 2p₋₁α|
 |1s₀β 2s₀α 2p₋₁β|
 |1s₀β 2s₀α 2p₀α|
 |1s₀β 2s₀α 2p₀β|
 |1s₀β 2s₀α 2p₁α|
 |1s₀β 2s₀α 2p₁β|
 |1s₀β 2s₀β 2p₋₁α|
 |1s₀β 2s₀β 2p₋₁β|
 |1s₀β 2s₀β 2p₀α|
 |1s₀β 2s₀β 2p₀β|
 |1s₀β 2s₀β 2p₁α|
 |1s₀β 2s₀β 2p₁β|
```

The operator we choose is the [Coulomb
repulsion](https://en.wikipedia.org/wiki/Coulomb%27s_law), implemented
by [`CoulombInteraction`](@ref):

```jldoctest lithium
julia> H = CoulombInteraction()
ĝ
```

- in case the configurations are identical, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ij}{ij}$:
  ```jldoctest lithium
  julia> Matrix(H, li)[1,1]
  F(1s₀α,1s₀β) - G(1s₀α,2s₀α) + F(1s₀α,2s₀α) + F(1s₀β,2s₀α)
  ```
  NB that some terms in the sum vanish due to spin-conservation in
  the two-body integral.
- in case the configurations differ by one orbital, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ik}{jk}$:
  ```jldoctest lithium
  julia> Matrix(H, [li[1],li_exc[6]])[2,1]
  - [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀α 2p₁β|2s₀α 1s₀β]
  ```
- in case the configurations differ by two orbitals, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ij}{kl}$:
  ```jldoctest lithium
  julia> Matrix(H, [li[1],li_exc[11]])[1,2]
  [1s₀β 2s₀α|2s₀β 2p₁α]
  ```
- in case the configurations differ by more than two orbital, $\twobodydx{\Phi_A}{\Phi_B} = 0$:
  ```jldoctest lithium
  julia> Matrix(H, [li[1],li_exc[23]])[1,2]
  0
  ```
  In this particular case, the matrix element vanishes because of
  spin-conservation, as well.

Again, we can generate the full two-body matrix:

```jldoctest lithium
julia> Matrix(H, vcat(li,li_exc))
26×26 Array{EnergyExpressions.NBodyMatrixElement,2}:
 F(1s₀α,1s₀β) - G(1s₀α,2s₀α) + F(1s₀α,2s₀α) + F(1s₀β,2s₀α)                 0                                                                           …  0
 0                                                                         F(1s₀α,1s₀β) + F(1s₀α,2s₀β) - G(1s₀β,2s₀β) + F(1s₀β,2s₀β)                      0
 0                                                                         0                                                                              0
 - [1s₀α 2p₋₁β|1s₀α 1s₀β] - [2s₀α 2p₋₁β|2s₀α 1s₀β]                         0                                                                              0
 0                                                                         0                                                                              0
 - [1s₀α 2p₀β|1s₀α 1s₀β] - [2s₀α 2p₀β|2s₀α 1s₀β]                           0                                                                           …  0
 0                                                                         0                                                                              0
 - [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀α 2p₁β|2s₀α 1s₀β]                           0                                                                              0
 [2s₀β 2p₋₁α|1s₀β 2s₀α]                                                    0                                                                              0
 0                                                                         - [1s₀α 2p₋₁β|1s₀α 1s₀β] - [2s₀β 2p₋₁β|2s₀β 1s₀β] + [2s₀β 2p₋₁β|1s₀β 2s₀β]     0
 [2s₀β 2p₀α|1s₀β 2s₀α]                                                     0                                                                           …  0
 0                                                                         - [1s₀α 2p₀β|1s₀α 1s₀β] - [2s₀β 2p₀β|2s₀β 1s₀β] + [2s₀β 2p₀β|1s₀β 2s₀β]        0
 [2s₀β 2p₁α|1s₀β 2s₀α]                                                     0                                                                              0
 0                                                                         - [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀β 2p₁β|2s₀β 1s₀β] + [2s₀β 2p₁β|1s₀β 2s₀β]        0
 [1s₀β 2p₋₁α|1s₀β 1s₀α] + [2s₀α 2p₋₁α|2s₀α 1s₀α] - [2s₀α 2p₋₁α|1s₀α 2s₀α]  0                                                                              0
 0                                                                         - [2s₀α 2p₋₁β|1s₀α 2s₀β]                                                    …  0
 [1s₀β 2p₀α|1s₀β 1s₀α] + [2s₀α 2p₀α|2s₀α 1s₀α] - [2s₀α 2p₀α|1s₀α 2s₀α]     0                                                                              0
 0                                                                         - [2s₀α 2p₀β|1s₀α 2s₀β]                                                        0
 [1s₀β 2p₁α|1s₀β 1s₀α] + [2s₀α 2p₁α|2s₀α 1s₀α] - [2s₀α 2p₁α|1s₀α 2s₀α]     0                                                                              0
 0                                                                         - [2s₀α 2p₁β|1s₀α 2s₀β]                                                        0
 0                                                                         [1s₀β 2p₋₁α|1s₀β 1s₀α] + [2s₀β 2p₋₁α|2s₀β 1s₀α]                             …  0
 0                                                                         0                                                                              - [1s₀β 2p₋₁β|2p₁β 1s₀β] + [1s₀β 2p₋₁β|1s₀β 2p₁β] - [2s₀β 2p₋₁β|2p₁β 2s₀β] + [2s₀β 2p₋₁β|2s₀β 2p₁β]
 0                                                                         [1s₀β 2p₀α|1s₀β 1s₀α] + [2s₀β 2p₀α|2s₀β 1s₀α]                                  0
 0                                                                         0                                                                              - [1s₀β 2p₀β|2p₁β 1s₀β] + [1s₀β 2p₀β|1s₀β 2p₁β] - [2s₀β 2p₀β|2p₁β 2s₀β] + [2s₀β 2p₀β|2s₀β 2p₁β]
 0                                                                         [1s₀β 2p₁α|1s₀β 1s₀α] + [2s₀β 2p₁α|2s₀β 1s₀α]                                  0
 0                                                                         0                                                                           …  - G(1s₀β,2s₀β) + F(1s₀β,2s₀β) - G(1s₀β,2p₁β) + F(1s₀β,2p₁β) - G(2s₀β,2p₁β) + F(2s₀β,2p₁β)
```

### Linear combination of N-body operators

Finally, we may form a linear combination of different many-body
operators:

```jldoctest helium
julia> H = OneBodyHamiltonian() + CoulombInteraction()
```

which allows to generate the full energy-expression matrix
simultaneously:

```jldoctest helium
julia> Matrix(H, vcat(he,he_exc[1:2]))
3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (1s₀α|1s₀α) + (1s₀β|1s₀β) + F(1s₀α,1s₀β)  0                                                            (1s₀β|2p₋₁β) + [1s₀α 1s₀β|1s₀α 2p₋₁β]
 0                                         (1s₀α|1s₀α) + (2p₋₁α|2p₋₁α) - G(1s₀α,2p₋₁α) + F(1s₀α,2p₋₁α)  0
 (2p₋₁β|1s₀β) + [1s₀α 2p₋₁β|1s₀α 1s₀β]     0                                                            (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β) + F(1s₀α,2p₋₁β)
```

## Non-orthogonal case, Löwdin rules

This case is more complex and is invoked by providing a list of
[`OrbitalOverlap`](@ref)s designating those pairs of orbitals which
are chosen to be _non-orthogonal_. Beware that the computation
complexity increases factorially with the amount of
non-orthogonalities! It is thus not a good idea to choose all orbitals
to non-orthogonal to one another.

For simplicity, we now consider a set of symbolic orbitals: `a,b,c`,
where specify that `b` and `c` are non-orthogonal:

```jldoctest non-orthogonal
julia> cfgs = SlaterDeterminant.([[:a, :b, :c], [:b, :d, :e]])
2-element Array{SlaterDeterminant{Symbol},1}:
 |a b c|
 |b d e|
```
The pairwise non-orthogonality can be specified simply as
```jldoctest non-orthogonal
julia> overlaps = [OrbitalOverlap(:b, :c)]
1-element Array{OrbitalOverlap{Symbol,Symbol},1}:
 ⟨b|c⟩
```

### One-body energy term

$$\onebody{\Phi_A}{\Phi_B} = (-)^{i+j}\onebody{i}{j}D^{AB}(i|j),$$

where $D^{AB}(i|j)$ is the determinant minor, where the row `i` and
column `j` are stricken out.

```jldoctest non-orthogonal
julia> Matrix(OneBodyHamiltonian(), cfgs, overlaps)
2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (a|a) - (a|a)⟨c|b⟩⟨b|c⟩ + (b|b) - (b|c)⟨c|b⟩ - (c|b)⟨b|c⟩ + (c|c)  0
 0                                                                  (b|b) + (d|d) + (e|e)
```

### Two-body energy term

Similarly, we have

$$\matrixel{\Phi_A}{\Omega_2}{\Phi_B} = (-)^{i+j+k+l}\matrixel{ij}{\Omega_2}{kl}D^{AB}(ij|kl),$$

where $D^{AB}(ij|kl)$ is the determinant minor, where the rows `i,j` and
columns `k,l` are stricken out.

The expressions become rather lengthy, so we only look at one
particular matrix element:

```jldoctest non-orthogonal
julia> Matrix(CoulombInteraction(), cfgs, overlaps)[1,2]
- ⟨c|b⟩[a b|e d] + ⟨c|b⟩[a b|d e] + [a c|e d] - [a c|d e]
```

### Higher-order terms

For details of the implementation and the general N-body case, see
[N-body matrix elements](@ref).

### References

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. I. Physical Interpretations by Means of Density Matrices,
  Natural Spin-Orbitals, and Convergence Problems in the Method of
  Configurational Interaction. Physical Review, 97(6),
  1474–1489. [10.1103/physrev.97.1474](http://dx.doi.org/10.1103/physrev.97.1474)

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. II. Study of the Ordinary Hartree-Fock
  Approximation. Physical Review, 97(6),
  1490–1508. [10.1103/physrev.97.1490](http://dx.doi.org/10.1103/physrev.97.1490)

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. III. Extension of the Hartree-Fock Scheme to Include
  Degenerate Systems and Correlation Effects. Physical Review, 97(6),
  1509–1520. [10.1103/physrev.97.1509](http://dx.doi.org/10.1103/physrev.97.1509)


```@meta
DocTestSetup = nothing
```
