# * Variation of orbital overlaps

"""
    diff(ab::OrbitalOverlap, o::O)

Vary the orbital overlap ⟨a|b⟩ with respect to |o⟩.
"""
Base.diff(ab::OrbitalOverlap, o::O) where O =
    o == ab.b ? NBodyEquation(Conjugate(ab.a),
                              IdentityOperator{1}()) : 0

"""
    diff(ab::OrbitalOverlap, o::Conjugate{O})

Vary the orbital overlap ⟨a|b⟩ with respect to ⟨o|.
"""
Base.diff(ab::OrbitalOverlap, o::Conjugate{O}) where O =
    conj(o) == ab.a ? NBodyEquation(ab.b,
                                    IdentityOperator{1}()) : 0

# * Variation of orbital matrix elements

"""
    diff(ome::OrbitalMatrixElement, o::O)

Vary the orbital overlap ⟨abc...|Ω|xyz...⟩ with respect to |o⟩.
"""
function Base.diff(ome::OrbitalMatrixElement{N}, o::O) where {N,O}
    i = findfirst(isequal(o), ome.b)
    isnothing(i) && return 0
    NBodyEquation(Conjugate(ome.a[i]), contract(ome, complement(N, i)...))
end

"""
    diff(ome::OrbitalMatrixElement, o::Conjugate{O})

Vary the orbital overlap ⟨abc...|Ω|xyz...⟩ with respect to ⟨o|.
"""
function Base.diff(ome::OrbitalMatrixElement{N}, o::Conjugate{O}) where {N,O}
    i = findfirst(isequal(conj(o)), ome.a)
    isnothing(i) && return 0
    NBodyEquation(ome.b[i], contract(ome, complement(N, i)...))
end

# * Variation of NBodyMatrixElements

"""
    diff(me::NBodyMatrixElement, o::O)

Vary the [`NBodyMatrixElement`](@ref) `me` with respect to the orbital `o`.
"""
function Base.diff(me::NBodyMatrixElement, o::O) where O
    equations = Vector{NBodyEquation}()
    for t in me.terms
        for (i,f) in enumerate(t.factors)
            v = diff(f, o)
            iszero(v) && continue
            v = v * NBodyTerm(vcat(t.factors[1:(i-1)]...,
                                   t.factors[(i+1):end]...),
                                   t.coeff)
            push!(equations, v)
        end
    end
    LinearCombinationEquation(equations)
end

# * Variation of energy expression matrices

"""
    diff(E::Matrix{NBodyMatrixElement}, o::O)

Vary the matrix of [`NBodyMatrixElement`](@ref)s with respect to the
orbital `o`.

# Examples

```jldoctest
julia> E = Matrix(OneBodyHamiltonian()+CoulombInteraction(),
                  SlaterDeterminant.([[:a, :b], [:c, :d]]))
2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:
 (a|a) + (b|b) - G(a,b) + F(a,b)  - [a b|d c] + [a b|c d]
 - [c d|b a] + [c d|a b]          (c|c) + (d|d) - G(c,d) + F(c,d)

julia> diff(E, :a)
2×2 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 2 stored entries:
  [1, 1]  =  ⟨a|ĥ + -⟨b|[a|b] + ⟨a|[b|b]
  [2, 1]  =  -⟨d|[c|b] + ⟨c|[d|b]

julia> diff(E, Conjugate(:b))
2×2 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 2 stored entries:
  [1, 1]  =  ĥ|b⟩ + -[a|b]|a⟩ + [a|a]|b⟩
  [1, 2]  =  -[a|d]|c⟩ + [a|c]|d⟩
```
"""
function Base.diff(E::EM, o::O) where {EM<:EnergyExpression,O}
    m,n = size(E)
    eqm = spzeros(LinearCombinationEquation, m, n)
    for i = 1:m
        for j = 1:n
            eq = diff(E[i,j], o)
            iszero(eq) && continue
            eqm[i,j] = eq
        end
    end
    eqm
end
