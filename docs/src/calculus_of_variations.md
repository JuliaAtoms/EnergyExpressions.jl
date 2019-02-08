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

``` # jldoctest helium
julia> he = spin_configurations(c"1s2")[1]
1s²

julia> a,b = he.orbitals
2-element Array{SpinOrbital{Orbital{Int64}},1}:
 1s₀α
 1s₀β
```

First we find the one- and two-body energy expressions:

``` # jldoctest helium
julia> h = OneBodyEnergyExpression(he,he)
I(1s₀α) + I(1s₀β)

julia> HC = TwoBodyEnergyExpression(he,he)
[1s₀β 1s₀α||1s₀β 1s₀α]
```

We can now vary these expressions with respect to the different
spin-orbitals; by convention, we vary the expressions with respect to
the _conjugate_ orbitals, to derive equations for the _unconjugated_
ones (this is important for complex orbitals, which is the case when
studying time-dependent problems):

``` # jldoctest helium
julia> h.integrals
2-element Array{OneBodyIntegral{SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}}},1}:
 I(1s₀α)
 I(1s₀β)

julia> diff.(h.integrals, Ref(conj(a)))
2-element Array{OneBodyHamiltonian,1}:
 ĥ1s₀α
 ĥ0
```

This result means that the variation of the first integral in the
one-body expression yields the one-body Hamiltonian $\hamiltonian$
acting on the orbital `1s₀α`, whereas the second integral, not
containing the first orbital, varies to yield zero.

If we instead vary with respect to the second orbital, _unconjugated_
this time, we get

``` # jldoctest helium
julia> diff.(h.integrals, Ref(b))
2-element Array{OneBodyHamiltonian,1}:
 ĥ0
 1s₀β†ĥ
```

Similarly, for the two-body energy expression:

``` # jldoctest helium
julia> diff.(HC.integrals, Ref(conj(a)))
1-element Array{DirectExchangePotentials{SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}}},1}:
 [1s₀β||1s₀β]1s₀α
```

```@meta
DocTestSetup = nothing
```
