# EnergyExpressions.jl

[![Documentation][docs-stable-img]][docs-stable-url]
[![Documentation (dev)][docs-dev-img]][docs-dev-url]
[![GitHub Actions CI][ci-gha-img]][ci-gha-url]
[![CodeCov][codecov-img]][codecov-url]

Small library for setting up energy expressions for electronic systems
that can be used to derive differential equations via calculus of
variations. The mathematics are not tied to atomic systems (which are
spherically symmetrical), and most of the code in this library could
be used to derive equations for e.g. molecules. However, the
`AbstractOrbital` and `SpinOrbital` types are taken from the
[AtomicLevels.jl](https://github.com/JuliaAtoms/AtomicLevels.jl)
library at present.

When dealing with configurations of spin-orbitals, where _all_ quantum
numbers are specified, it is trivial to derive the energy
expression. For e.g. configuration state functions (CSFs), which are
spin-averaged linear combinations of multiple Slater determinants, the
problem is much more involved and requires dealing with angular
momentum coupling and tensor algebra. This library could be used as a
basis for building such functionality.

[ci-gha-url]: https://github.com/JuliaAtoms/EnergyExpressions.jl/actions
[ci-gha-img]: https://github.com/JuliaAtoms/EnergyExpressions.jl/workflows/CI/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaAtoms/EnergyExpressions.jl
[codecov-img]: https://codecov.io/gh/JuliaAtoms/EnergyExpressions.jl/branch/master/graph/badge.svg
[docs-stable-url]: https://juliaatoms.org/EnergyExpressions.jl/stable/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-url]: https://juliaatoms.org/EnergyExpressions.jl/dev/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
