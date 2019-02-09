# Calculus of Variations

[Calculus of
Variations](https://en.wikipedia.org/wiki/Calculus_of_variations) is a
mathematical technique for deriving differential equations, whose
solutions fulfil a certain criterion. For the case at hand, we are
interested in solutions whose total energy matches with the one set up
by the energy expression for all the electrons in a system. The
following variational rules are useful:

$$\begin{equation}
\vary{\bra{a}}\braket{a}{b} = \ket{b},
\end{equation}$$
where the notation $\vary{\bra{a}}\braket{a}{b}$ means _vary
$\braket{a}{b}$ with respect to $\bra{a}$_.

$$\begin{equation}
\begin{cases}
\vary{\bra{l}}\onebody{l}{k} &=
\vary{\bra{l}}\matrixel{l}{\hamiltonian}{k} =
\hamiltonian\ket{k},\\
\vary{\ket{k}}\onebody{l}{k} &=
\bra{l}\hamiltonian.
\end{cases}
\end{equation}$$

$$\begin{equation}
\begin{aligned}
\vary{\bra{a}}\twobodydx{ab}{cd} &=
\vary{\bra{a}}\{\twobody{ab}{cd}-\twobody{ab}{dc}\} =
\vary{\bra{a}}\{\matrixel{ab}{r_{12}^{-1}}{cd}-\matrixel{ab}{r_{12}^{-1}}{dc}\} \\
&=
\matrixel{b}{r_{12}^{-1}}{d}\ket{c}-
\matrixel{b}{r_{12}^{-1}}{c}\ket{d} \equiv
\twobody{b}{d}\ket{c}-
\twobody{b}{c}\ket{d}\defd
(\direct{bd}-\exchange{bd})\ket{c}.
\end{aligned}
\end{equation}$$

## Helium

Returning to our helium example (now only considering the ground
state):

```@meta
DocTestSetup = quote
    using AtomicLevels
    using EnergyExpressions
end
```

```jldoctest helium-variation
julia> he = SlaterDeterminant.(spin_configurations(c"1s2"))
1-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:
 |1s₀α 1s₀β|

julia> a,b = he[1].orbitals
2-element Array{SpinOrbital{Orbital{Int64}},1}:
 1s₀α
 1s₀β
```

First we find the energy expression:

```jldoctest helium-variation
julia> E = Matrix(OneBodyHamiltonian() + CoulombInteraction(), he)
1×1 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (1s₀α|1s₀α) + (1s₀β|1s₀β) + F(1s₀α,1s₀β)
```

We can now vary this expression with respect to the different
spin-orbitals; by convention, we vary the expression with respect to
the _conjugate_ orbitals, to derive equations for the _unconjugated_
ones (this is important for complex orbitals, which is the case when
studying time-dependent problems):

```jldoctest helium-variation
julia> diff(E, Conjugate(a))
1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:
  [1, 1]  =  ĥ|1s₀α⟩ + [1s₀β|1s₀β]|1s₀α⟩

julia> diff(E, Conjugate(b))
1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:
  [1, 1]  =  ĥ|1s₀β⟩ + [1s₀α|1s₀α]|1s₀β⟩
```

This result means that the variation of the first integral in the
one-body part of the energy expression yields the one-body Hamiltonian
$\hamiltonian$ acting on the orbital `1s₀α`, whereas the second
integral, not containing the first orbital, varies to yield zero. The
two-body energy `F(1s₀α,1s₀β)` varies to yield the two-body potential
`[1s₀β|1s₀β]`, again acting on the orbital `1s₀α`.

If we instead vary with respect to the _unconjugated_ orbitals, we get

```jldoctest helium-variation
julia> diff(E, a)
1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:
  [1, 1]  =  ⟨1s₀α|ĥ + ⟨1s₀α|[1s₀β|1s₀β]

julia> diff(E, b)
1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:
  [1, 1]  =  ⟨1s₀β|ĥ + ⟨1s₀β|[1s₀α|1s₀α]
```

that is, we instead derived equations for the conjugated orbitals.

```@meta
DocTestSetup = nothing
```
