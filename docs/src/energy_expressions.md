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

The examples below are for atomic configurations, but the library is
not fundamentally limited to this case.

## Orthogonal case, Slater–Condon rules

```@meta
DocTestSetup = quote
    using AtomicLevels
    using EnergyExpressions
end
```

First, we find all configurations possible for two configurations of
helium:

```jldoctest helium
julia> he = spin_configurations(c"1s2")
1-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:
 1s²

julia> he_exc = spin_configurations(c"1s 2p")
12-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:
 1s₀α 2p₋₁α
 1s₀α 2p₋₁β
 1s₀α 2p₀α
 1s₀α 2p₀β
 1s₀α 2p₁α
 1s₀α 2p₁β
 1s₀β 2p₋₁α
 1s₀β 2p₋₁β
 1s₀β 2p₀α
 1s₀β 2p₀β
 1s₀β 2p₁α
 1s₀β 2p₁β
```

### One-body energy term

The Slater–Condon rules state that the one-body energy expression between two
configurations is

- in case the configurations are identical, $\onebody{\Phi_A}{\Phi_B} = \onebody{i}{i}$:
  ```jldoctest helium
  julia> OneBodyEnergyExpression(he[1],he[1])
  I(1s₀α) + I(1s₀β)
  ```
- in case the configurations differ by one orbital, $\onebody{\Phi_A}{\Phi_B} = \onebody{i}{j}$:
  ```jldoctest helium
  julia> OneBodyEnergyExpression(he[1],he_exc[1])
  I(1s₀β, 2p₋₁α)
  ```
- in case the configurations differ by more than one orbital, $\onebody{\Phi_A}{\Phi_B} = 0$:
  ```jldoctest helium
  julia> OneBodyEnergyExpression(he_exc[1],he_exc[9])
  0
  ```

We can easily generate the full one-body matrix:

```jldoctest helium
julia> one_body_hamiltonian_matrix(SpinOrbital, he_exc)
12×12 Array{OneBodyEnergyExpression{SpinOrbital,SpinOrbital},2}:
 I(1s₀α) + I(2p₋₁α)  I(2p₋₁α, 2p₋₁β)     I(2p₋₁α, 2p₀α)     I(2p₋₁α, 2p₀β)     I(2p₋₁α, 2p₁α)     I(2p₋₁α, 2p₁β)     I(1s₀α, 1s₀β)       0                   0                  0                  0                  0
 I(2p₋₁β, 2p₋₁α)     I(1s₀α) + I(2p₋₁β)  I(2p₋₁β, 2p₀α)     I(2p₋₁β, 2p₀β)     I(2p₋₁β, 2p₁α)     I(2p₋₁β, 2p₁β)     0                   I(1s₀α, 1s₀β)       0                  0                  0                  0
 I(2p₀α, 2p₋₁α)      I(2p₀α, 2p₋₁β)      I(1s₀α) + I(2p₀α)  I(2p₀α, 2p₀β)      I(2p₀α, 2p₁α)      I(2p₀α, 2p₁β)      0                   0                   I(1s₀α, 1s₀β)      0                  0                  0
 I(2p₀β, 2p₋₁α)      I(2p₀β, 2p₋₁β)      I(2p₀β, 2p₀α)      I(1s₀α) + I(2p₀β)  I(2p₀β, 2p₁α)      I(2p₀β, 2p₁β)      0                   0                   0                  I(1s₀α, 1s₀β)      0                  0
 I(2p₁α, 2p₋₁α)      I(2p₁α, 2p₋₁β)      I(2p₁α, 2p₀α)      I(2p₁α, 2p₀β)      I(1s₀α) + I(2p₁α)  I(2p₁α, 2p₁β)      0                   0                   0                  0                  I(1s₀α, 1s₀β)      0
 I(2p₁β, 2p₋₁α)      I(2p₁β, 2p₋₁β)      I(2p₁β, 2p₀α)      I(2p₁β, 2p₀β)      I(2p₁β, 2p₁α)      I(1s₀α) + I(2p₁β)  0                   0                   0                  0                  0                  I(1s₀α, 1s₀β)
 I(1s₀β, 1s₀α)       0                   0                  0                  0                  0                  I(1s₀β) + I(2p₋₁α)  I(2p₋₁α, 2p₋₁β)     I(2p₋₁α, 2p₀α)     I(2p₋₁α, 2p₀β)     I(2p₋₁α, 2p₁α)     I(2p₋₁α, 2p₁β)
 0                   I(1s₀β, 1s₀α)       0                  0                  0                  0                  I(2p₋₁β, 2p₋₁α)     I(1s₀β) + I(2p₋₁β)  I(2p₋₁β, 2p₀α)     I(2p₋₁β, 2p₀β)     I(2p₋₁β, 2p₁α)     I(2p₋₁β, 2p₁β)
 0                   0                   I(1s₀β, 1s₀α)      0                  0                  0                  I(2p₀α, 2p₋₁α)      I(2p₀α, 2p₋₁β)      I(1s₀β) + I(2p₀α)  I(2p₀α, 2p₀β)      I(2p₀α, 2p₁α)      I(2p₀α, 2p₁β)
 0                   0                   0                  I(1s₀β, 1s₀α)      0                  0                  I(2p₀β, 2p₋₁α)      I(2p₀β, 2p₋₁β)      I(2p₀β, 2p₀α)      I(1s₀β) + I(2p₀β)  I(2p₀β, 2p₁α)      I(2p₀β, 2p₁β)
 0                   0                   0                  0                  I(1s₀β, 1s₀α)      0                  I(2p₁α, 2p₋₁α)      I(2p₁α, 2p₋₁β)      I(2p₁α, 2p₀α)      I(2p₁α, 2p₀β)      I(1s₀β) + I(2p₁α)  I(2p₁α, 2p₁β)
 0                   0                   0                  0                  0                  I(1s₀β, 1s₀α)      I(2p₁β, 2p₋₁α)      I(2p₁β, 2p₋₁β)      I(2p₁β, 2p₀α)      I(2p₁β, 2p₀β)      I(2p₁β, 2p₁α)      I(1s₀β) + I(2p₁β)
```

### Two-body energy term

Similar considerations apply for the two-body energy terms between two
configurations. To make it more interesting, we consider lithium
which has three electrons:

```jldoctest lithium
julia> li = spin_configurations(c"1s2 2s")
2-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:
 1s² 2s₀α
 1s² 2s₀β

julia> li_exc = spin_configurations(c"1s 2s 2p")
24-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:
 1s₀α 2s₀α 2p₋₁α
 1s₀α 2s₀α 2p₋₁β
 1s₀α 2s₀α 2p₀α
 1s₀α 2s₀α 2p₀β
 1s₀α 2s₀α 2p₁α
 1s₀α 2s₀α 2p₁β
 1s₀α 2s₀β 2p₋₁α
 1s₀α 2s₀β 2p₋₁β
 1s₀α 2s₀β 2p₀α
 1s₀α 2s₀β 2p₀β
 1s₀α 2s₀β 2p₁α
 1s₀α 2s₀β 2p₁β
 1s₀β 2s₀α 2p₋₁α
 1s₀β 2s₀α 2p₋₁β
 1s₀β 2s₀α 2p₀α
 1s₀β 2s₀α 2p₀β
 1s₀β 2s₀α 2p₁α
 1s₀β 2s₀α 2p₁β
 1s₀β 2s₀β 2p₋₁α
 1s₀β 2s₀β 2p₋₁β
 1s₀β 2s₀β 2p₀α
 1s₀β 2s₀β 2p₀β
 1s₀β 2s₀β 2p₁α
 1s₀β 2s₀β 2p₁β
```

- in case the configurations are identical, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ij}{ij}$:
  ```jldoctest lithium
  julia> TwoBodyEnergyExpression(li[1],li[1])
  [1s₀β 1s₀α||1s₀β 1s₀α] + [2s₀α 1s₀α||2s₀α 1s₀α] + [2s₀α 1s₀β||2s₀α 1s₀β]
  ```
- in case the configurations differ by one orbital, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ik}{jk}$:
  ```jldoctest lithium
  julia> TwoBodyEnergyExpression(li[1],li[2])
  [2s₀α 1s₀α||2s₀β 1s₀α] + [2s₀α 1s₀β||2s₀β 1s₀β]
  ```
- in case the configurations differ by two orbitals, $\twobodydx{\Phi_A}{\Phi_B} = \twobodydx{ij}{kl}$:
  ```jldoctest lithium
  julia> TwoBodyEnergyExpression(li[1],li_exc[end])
  [1s₀α 2s₀α||2s₀β 2p₁β]
  ```
  
- in case the configurations differ by more than two orbital, $\twobodydx{\Phi_A}{\Phi_B} = 0$:
  ```jldoctest lithium
  julia> TwoBodyEnergyExpression(li_exc[1],li_exc[end])
  0
  ```

Again, we can generate the full two-body matrix:

```jldoctest lithium
julia> two_body_hamiltonian_matrix(SpinOrbital, li_exc)
24×24 Array{TwoBodyEnergyExpression{SpinOrbital,SpinOrbital,SpinOrbital,SpinOrbital},2}:
 [2s₀α 1s₀α||2s₀α 1s₀α] + … + [2p₋₁α 2s₀α||2p₋₁α 2s₀α]  [2p₋₁α 1s₀α||2p₋₁β 1s₀α] + [2p₋₁α 2s₀α||2p₋₁β 2s₀α]    …  0
 [2p₋₁β 1s₀α||2p₋₁α 1s₀α] + [2p₋₁β 2s₀α||2p₋₁α 2s₀α]    [2s₀α 1s₀α||2s₀α 1s₀α] + … + [2p₋₁β 2s₀α||2p₋₁β 2s₀α]     0
 [2p₀α 1s₀α||2p₋₁α 1s₀α] + [2p₀α 2s₀α||2p₋₁α 2s₀α]      [2p₀α 1s₀α||2p₋₁β 1s₀α] + [2p₀α 2s₀α||2p₋₁β 2s₀α]         0
 [2p₀β 1s₀α||2p₋₁α 1s₀α] + [2p₀β 2s₀α||2p₋₁α 2s₀α]      [2p₀β 1s₀α||2p₋₁β 1s₀α] + [2p₀β 2s₀α||2p₋₁β 2s₀α]         0
 [2p₁α 1s₀α||2p₋₁α 1s₀α] + [2p₁α 2s₀α||2p₋₁α 2s₀α]      [2p₁α 1s₀α||2p₋₁β 1s₀α] + [2p₁α 2s₀α||2p₋₁β 2s₀α]         0
 [2p₁β 1s₀α||2p₋₁α 1s₀α] + [2p₁β 2s₀α||2p₋₁α 2s₀α]      [2p₁β 1s₀α||2p₋₁β 1s₀α] + [2p₁β 2s₀α||2p₋₁β 2s₀α]      …  [1s₀α 2s₀α||1s₀β 2s₀β]
 [2s₀β 1s₀α||2s₀α 1s₀α] + [2s₀β 2p₋₁α||2s₀α 2p₋₁α]      [2s₀β 2p₋₁α||2s₀α 2p₋₁β]                                  [1s₀α 2p₋₁α||1s₀β 2p₁β]
 [2s₀β 2p₋₁β||2s₀α 2p₋₁α]                               [2s₀β 1s₀α||2s₀α 1s₀α] + [2s₀β 2p₋₁β||2s₀α 2p₋₁β]         [1s₀α 2p₋₁β||1s₀β 2p₁β]
 [2s₀β 2p₀α||2s₀α 2p₋₁α]                                [2s₀β 2p₀α||2s₀α 2p₋₁β]                                   [1s₀α 2p₀α||1s₀β 2p₁β]
 [2s₀β 2p₀β||2s₀α 2p₋₁α]                                [2s₀β 2p₀β||2s₀α 2p₋₁β]                                   [1s₀α 2p₀β||1s₀β 2p₁β]
 [2s₀β 2p₁α||2s₀α 2p₋₁α]                                [2s₀β 2p₁α||2s₀α 2p₋₁β]                                …  [1s₀α 2p₁α||1s₀β 2p₁β]
 [2s₀β 2p₁β||2s₀α 2p₋₁α]                                [2s₀β 2p₁β||2s₀α 2p₋₁β]                                   [1s₀α 2s₀β||1s₀β 2s₀β] + [1s₀α 2p₁β||1s₀β 2p₁β]
 [1s₀β 2s₀α||1s₀α 2s₀α] + [1s₀β 2p₋₁α||1s₀α 2p₋₁α]      [1s₀β 2p₋₁α||1s₀α 2p₋₁β]                                  [2s₀α 2p₋₁α||2s₀β 2p₁β]
 [1s₀β 2p₋₁β||1s₀α 2p₋₁α]                               [1s₀β 2s₀α||1s₀α 2s₀α] + [1s₀β 2p₋₁β||1s₀α 2p₋₁β]         [2s₀α 2p₋₁β||2s₀β 2p₁β]
 [1s₀β 2p₀α||1s₀α 2p₋₁α]                                [1s₀β 2p₀α||1s₀α 2p₋₁β]                                   [2s₀α 2p₀α||2s₀β 2p₁β]
 [1s₀β 2p₀β||1s₀α 2p₋₁α]                                [1s₀β 2p₀β||1s₀α 2p₋₁β]                                …  [2s₀α 2p₀β||2s₀β 2p₁β]
 [1s₀β 2p₁α||1s₀α 2p₋₁α]                                [1s₀β 2p₁α||1s₀α 2p₋₁β]                                   [2s₀α 2p₁α||2s₀β 2p₁β]
 [1s₀β 2p₁β||1s₀α 2p₋₁α]                                [1s₀β 2p₁β||1s₀α 2p₋₁β]                                   [2s₀α 1s₀β||2s₀β 1s₀β] + [2s₀α 2p₁β||2s₀β 2p₁β]
 [1s₀β 2s₀β||1s₀α 2s₀α]                                 0                                                         [2p₋₁α 1s₀β||2p₁β 1s₀β] + [2p₋₁α 2s₀β||2p₁β 2s₀β]
 0                                                      [1s₀β 2s₀β||1s₀α 2s₀α]                                    [2p₋₁β 1s₀β||2p₁β 1s₀β] + [2p₋₁β 2s₀β||2p₁β 2s₀β]
 0                                                      0                                                      …  [2p₀α 1s₀β||2p₁β 1s₀β] + [2p₀α 2s₀β||2p₁β 2s₀β]
 0                                                      0                                                         [2p₀β 1s₀β||2p₁β 1s₀β] + [2p₀β 2s₀β||2p₁β 2s₀β]
 0                                                      0                                                         [2p₁α 1s₀β||2p₁β 1s₀β] + [2p₁α 2s₀β||2p₁β 2s₀β]
 0                                                      0
 [2s₀β 1s₀β||2s₀β 1s₀β] + … + [2p₁β 2s₀β||2p₁β 2s₀β]
```

## Non-orthogonal case, Löwdin rules

This case is more complex and not yet fully implemented in
EnergyExpressions.jl

### One-body energy term

$$\onebody{\Phi_A}{\Phi_B} = (-)^{i+j}\onebody{i}{j}D^{AB}_{ji},$$

where $D^{AB}_{ji}$ is the determinant minor.

The implementation is only correct for two-electron configurations at
the moment:

```jldoctest helium
julia> OneBodyEnergyExpression{SpinOrbital,SpinOrbital}(he[1],he[1],orthogonal=false)
I(1s₀α) - I(1s₀α, 1s₀β) - I(1s₀β, 1s₀α) + I(1s₀β)
```

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
