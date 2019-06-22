# * N-body term factors

"""
    NBodyTermFactor

Abstract type for a factor in a term in a N-body matrix element expansion
"""
abstract type NBodyTermFactor end

# ** Orbital overlaps

"""
    OrbitalOverlap(a,b)

Represents the overlap between the orbitals `a` and `b` in a N-body matrix element expansion.

# Examples

```jldoctest
julia> EnergyExpressions.OrbitalOverlap(:a,:b)
⟨a|b⟩
```
"""
struct OrbitalOverlap{A,B} <: NBodyTermFactor
    a::A
    b::B
end

# By default, we assume that all orbital overlaps are non-zero,
# i.e. the orbitals are non-orthogonal, since we specify
# non-orthogonality precisely by providing OrbitalOverlaps.
Base.iszero(::OrbitalOverlap) = false

Base.adjoint(o::OrbitalOverlap{A,B}) where {A,B} = OrbitalOverlap{B,A}(o.b, o.a)

"""
    numbodies(::OrbitalOverlap)

Returns the number of bodies coupled by the zero-body operator in the
orbital overlap, i.e. `0`.
"""
numbodies(::OrbitalOverlap) = 0

"""
    isdependent(o::OrbitalOverlap, orbital)

Returns `true` if the [`OrbitalOverlap`](@ref) `o` depends on `orbital`.

# Examples

```jldoctest
julia> isdependent(OrbitalOverlap(:a,:b), :a)
false

julia> isdependent(OrbitalOverlap(:a,:b), Conjugate(:a))
true

julia> isdependent(OrbitalOverlap(:a,:b), :b)
true
```
"""
isdependent(o::OrbitalOverlap, corb::Conjugate{O}) where O = o.a == corb.orbital
isdependent(o::OrbitalOverlap, orb::O) where O = o.b == orb

Base.show(io::IO, o::OrbitalOverlap) =
    write(io, "⟨$(o.a)|$(o.b)⟩")

# ** Orbital matrix elements

"""
    OrbitalMatrixElement(a,o,b)

Represents the N-body matrix element between the sets of orbitals `a` and `b`.

# Examples

```jldoctest
julia> struct MyTwoBodyOperator <: TwoBodyOperator end

julia> EnergyExpressions.OrbitalMatrixElement((:a,:b), MyTwoBodyOperator(), (:c,:d))
⟨a b|MyTwoBodyOperator()|c d⟩
```
"""
struct OrbitalMatrixElement{N,A,O<:NBodyOperator{N},B} <: NBodyTermFactor
    a::Vector{A}
    o::O
    b::Vector{B}
end

Base.iszero(::OrbitalMatrixElement) = false

Base.:(==)(a::OrbitalMatrixElement, b::OrbitalMatrixElement) =
    a.a == b.a && a.o == b.o && a.b == b.b

"""
    numbodies(::OrbitalMatrixElement{N})

Returns the number of bodies coupled by the operator, i.e. `N`.
"""
numbodies(::OrbitalMatrixElement{N}) where N = N

"""
    contract(ome::OrbitalMatrixElement{N}, i...)

Contract `ome` over all coordinates `i...`. `length(i)` cannot be
larger than `N`.
"""
function contract(ome::OrbitalMatrixElement{N,A,O,B}, i::Integer...) where {N,A,O,B}
    Q = N - length(i)

    Q > 0 ||
        throw(ArgumentError("Cannot contract $(N)-body orbital matrix element over $(length(i)) coordinates"))

    Q == N && return ome.o

    ContractedOperator(NTuple{Q,A}(ome.a[ii] for ii in i),
                       ome.o,
                       NTuple{Q,B}(ome.b[ii] for ii in i))
end


"""
    isdependent(o::OrbitalMatrixElement, orbital)

Returns `true` if the [`OrbitalMatrixElement`](@ref) `o` depends on `orbital`.

# Examples

```jldoctest
julia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), :a)
false

julia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), Conjugate(:a))
true

julia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), :b)
true

julia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,:b,), CoulombInteraction(), (:c,:d)), :c)
true
```
"""
isdependent(o::OrbitalMatrixElement, corb::Conjugate{O}) where O = corb.orbital ∈ o.a
isdependent(o::OrbitalMatrixElement, orb::O) where O = orb ∈ o.b

function Base.show(io::IO, ome::OrbitalMatrixElement{N}) where N
    write(io, "⟨", join(string.(ome.a), " "))
    N > 0 && write(io, "|")
    show(io, ome.o)
    N > 0 && write(io, "|")
    write(io, join(string.(ome.b), " "), "⟩")
end

# * N-body terms

"""
    NBodyTerm(factors, coeff)

Structure representing one term in the expansion of a N-body matrix element.
"""
struct NBodyTerm
    factors::Vector{NBodyTermFactor}
    coeff
end

Base.one(::Type{NBodyTerm}) = NBodyTerm(NBodyTermFactor[], 1)
Base.one(::NBodyTerm) = one(NBodyTerm)
Base.isone(term::NBodyTerm) = isempty(term.factors) && isone(term.coeff)
isminusone(term::NBodyTerm) = isempty(term.factors) && isone(-term.coeff)

Base.zero(::Type{NBodyTerm}) = NBodyTerm(NBodyTermFactor[], 0)
Base.zero(::NBodyTerm) = zero(NBodyTerm)
Base.iszero(term::NBodyTerm) = iszero(term.coeff) || any(iszero.(term.factors))

Base.convert(::Type{NBodyTerm}, f::NBodyTermFactor) =
    NBodyTerm([f], 1)

Base.:(*)(a::NBodyTerm, b::NBodyTerm) =
    NBodyTerm(vcat(a.factors, b.factors), a.coeff*b.coeff)

Base.:(*)(a::NBodyTerm, b::Number) =
    NBodyTerm(a.factors, a.coeff*b)

Base.:(*)(a::Number, b::NBodyTerm) =
    NBodyTerm(b.factors, a*b.coeff)

Base.:(*)(a::NBodyTerm, b::NBodyTermFactor) =
    NBodyTerm(vcat(a.factors, b), a.coeff)

Base.:(*)(a::NBodyTermFactor, b::NBodyTerm) =
    NBodyTerm(vcat(a, b.factors), b.coeff)

Base.:(*)(a::NBodyTermFactor, b::NBodyTermFactor) =
    NBodyTerm([a, b], 1)

Base.:(*)(a::NBodyTermFactor, b::Number) =
    NBodyTerm([a], b)

Base.:(*)(a::Number, b::NBodyTermFactor) =
    NBodyTerm([b], a)

"""
    isdependent(nbt::NBodyTerm, o)

Returns `true` if any of the factors comprising `nbt` is dependent on
the orbital `o`. Not that the result is dependent on whether `o` is
conjugated or not.

"""
isdependent(nbt::NBodyTerm, o) =
    any(isdependent.(nbt.factors, Ref(o)))

"""
    transform(f::Function, nbt::NBodyTerm)

Transform integrals of the the N-body matrix element expansion term
`nbt` according to the function `f`, which should accept a single
[`NBodyTermFactor`](@ref) as its argument.
"""
function transform(f::Function, nbt::NBodyTerm)
    isempty(nbt.factors) && return NBodyMatrixElement([nbt])

    # Generate the starting N-body matrix element by applying the
    # transform function `f` on the first factor of the NBodyTerm nbt.
    nbme = [nbt.coeff*f(first(nbt.factors))]

    # Then generate successively larger expansion by multiplying all
    # previous terms by the transformation of the current factor.
    for fac in nbt.factors[2:end]
        nbme[1] = sum([term*f(fac) for term in nbme[1].terms])
    end

    nbme[1]
end

function Base.show(io::IO, term::NBodyTerm; show_sign=false)
    showcoeff(io, term.coeff, show_sign, isempty(term.factors))
    (iszero(term) || isempty(term.factors)) && return
    noprintfactors = 2
    factors = if get(io, :limit, false) && length(term.factors) > 2noprintfactors
        vcat(term.factors[1:noprintfactors], "…", term.factors[end-(noprintfactors-1):end])
    else
        term.factors
    end
    write(io, join(string.(factors), ""))
end

# * N-body matrix elements

"""
    NBodyMatrixElement(terms)

Structure representing the expansion of a N-body matrix element.
"""
struct NBodyMatrixElement
    terms::Vector{NBodyTerm}
end

Base.zero(::Type{NBodyMatrixElement}) =
    NBodyMatrixElement(NBodyTerm[])

Base.zero(::NBodyMatrixElement) = zero(NBodyMatrixElement)

Base.one(::Type{NBodyMatrixElement}) =
    NBodyMatrixElement([one(NBodyTerm)])

Base.iszero(nbme::NBodyMatrixElement) =
    isempty(nbme.terms) || all(iszero.(nbme.terms))

Base.convert(::Type{NBodyMatrixElement}, nbt::NBodyTerm) =
    NBodyMatrixElement([nbt])

Base.:(+)(a::NBodyMatrixElement, b::NBodyMatrixElement) =
    NBodyMatrixElement(vcat(a.terms, b.terms))

Base.:(+)(a::NBodyMatrixElement, b::NBodyTerm) =
    NBodyMatrixElement(vcat(a.terms, b))

Base.:(+)(a::NBodyTerm, b::NBodyMatrixElement) =
    NBodyMatrixElement(vcat(a, b.terms))

Base.:(+)(a::NBodyTerm, b::NBodyTerm) =
    NBodyMatrixElement([a, b])

Base.:(*)(a::NBodyMatrixElement, b::Union{Number,NBodyTerm,NBodyTermFactor}) =
    NBodyMatrixElement([b*t for t in a.terms])

Base.:(*)(a::Union{Number,NBodyTerm,NBodyTermFactor}, b::NBodyMatrixElement) =
    NBodyMatrixElement([a*t for t in b.terms])

function Base.show(io::IO, nbme::NBodyMatrixElement)
    if iszero(nbme)
        write(io, "0")
        return
    end
    noprintterms = 6
    nterms = length(nbme.terms)
    for (i,term) in enumerate(nbme.terms)
        if get(io, :limit, false) && i > noprintterms && i ≤ nterms-noprintterms
            i == noprintterms+1 && write(io, " + …")
        else
            i > 1 && write(io, " ")
            show(io, term, show_sign=(i > 1))
        end
    end
end

# ** Determinant utilities

"""
    detaxis(i::CartesianIndex{N})

Generate the axis index vector for the determinant minor, whose rows
or columns represented by the `CartesianIndex` `i` should be omitted.
Implemented via [`complement`](@ref).
"""
detaxis(i::AbstractVector, n) = complement(n, i)
detaxis(i::Tuple, n) = complement(n, i...)
detaxis(i::CartesianIndex, n) = complement(n, Tuple(i)...)
detaxis(i::Integer, n) = complement(n, i)

"""
    detminor(k, l, A)

Calculate the [determinant
minor](https://en.wikipedia.org/wiki/Minor_(linear_algebra)) of `A`,
where the rows `k` and the columns `l` have been stricken out.
"""
function detminor(k, l, A)
    n = size(A,1)
    n == 1 && return 1
    det(A[detaxis(k,n),detaxis(l,n)])
end

indexsum(k::Integer,l::Integer) = k+l
indexsum(k::Tuple,l::Tuple) = sum(k) + sum(l)
indexsum(k::AbstractVector,l::AbstractVector) = sum(k) + sum(l)
indexsum(k::CartesianIndex,l::CartesianIndex) = indexsum(Tuple(k), Tuple(l))

"""
    powneg1(k) = (-)ᵏ

Calculates powers of negative unity for integer `k`.
"""
powneg1(k::Integer) = isodd(k) ? -1 : 1

"""
    cofactor(k, l, A)

Calculate the
[cofactor](https://en.wikipedia.org/wiki/Minor_(linear_algebra)) of
`A`, where the rows `k` and the columns `l` have been stricken out.
The cofactor is calculated recursively, by expanding the minor
determinants in cofactors, so this function should only be used in
case it is known that the cofactor is non-zero.
"""
cofactor(k, l, A) = powneg1(indexsum(k,l))*detminor(k, l, A)

function merge_overlapping_blocks(blocks::Vector{UnitRange{Int}})
    cur_block = first(blocks)
    new_blocks = UnitRange{Int}[]
    for b in blocks
        if isempty(cur_block ∩ b)
            push!(new_blocks, cur_block)
            cur_block = b
        else
            cur_block = min(cur_block[1],b[1]):max(cur_block[end],b[end])
        end
    end
    if isempty(new_blocks)
        push!(new_blocks, cur_block)
    elseif blocks[end] ∉ new_blocks[end]
        push!(new_blocks, blocks[end])
    end

    new_blocks
end

function find_diagonal_blocks(A::AbstractSparseMatrix)
    m,n = size(A)
    m == n || throw(ArgumentError("Matrix must be square"))

    rows = rowvals(A)

    # Find the row extent of each column, including the diagonal; each
    # column is tentatively its own diagonal block.
    blocks = map(1:n) do j
        j_rows = rows[nzrange(A,j)]
        min(j,minimum(j_rows)):max(j,maximum(j_rows))
    end

    # Merge overlapping blocks; since they are ordered by column and
    # always contain the diagonal, we can do a forward sweep and then
    # a backward sweep to capture all.
    blocks = merge_overlapping_blocks(blocks)
    reverse(merge_overlapping_blocks(reverse(blocks)))
end

function block_det(A::AbstractSparseMatrix{T}, blocks::Vector{UnitRange{Int}}) where T
    D = one(T)
    for b in blocks
        if length(b) > 1
            D *= det(A[b,b])
        else
            i = b[1]
            D *= A[i,i]
        end
    end
    D
end

"""
    det(A)

Calculate the determinant of the matrix `A` whose elements are of the
[`NBodyTerm`](@ref) type, by expanding the determinant along the first
column. This is an expensive operation, and should only be done with
relatively sparse matrices.
"""
function LinearAlgebra.det(A::AbstractSparseMatrix{<:NBodyTerm})
    m,n = size(A)
    # Some trivial special cases
    (m,n) == (0,0) && return 1
    (m,n) == (1,1) && return A[1,1]
    iszero(A) && return zero(NBodyMatrixElement)

    # If the matrix is too large, we try to break it up in smaller
    # diagonal blocks; however, if we only get one block back, we have
    # to apply the general case.
    if m > 5
        blocks = find_diagonal_blocks(A)
        length(blocks) > 1 && return block_det(A, blocks)
    end

    # The general case, Laplace expansion.
    D = zero(NBodyMatrixElement)
    j = 1 # Expand along first column
    for i = 1:m
        iszero(A[i,j]) && continue
        Cᵢⱼ = cofactor(i,j,A)
        iszero(Cᵢⱼ) && continue
        D += A[i,j]*Cᵢⱼ
    end
    D
end

"""
    permutation_sign(p)

Calculate the sign of the permutation `p`, 1 if `iseven(p)`, -1
otherwise.
"""
permutation_sign(p) =
    iseven(Combinatorics.parity(p)) ? 1 : -1

# ** N-body matrix construction

"""
    NBodyMatrixElement(a, op, b, overlap)

Generate the matrix element of the N-body operator `op`, between the
Slater determinants `a` and `b`, according to the Löwdin rules. The
matrix `overlap` contains the mutual overlaps between all
single-particle orbitals in the Slater determinants. If the orbitals
are all orthogonal, the Löwdin rules collapse to the Slater–Condon
rules.
"""
@generated function NBodyMatrixElement(a::SlaterDeterminant, op::NBodyOperator{N}, b::SlaterDeterminant, overlap) where N
    quote
        length(a) == length(b) ||
            throw(DimensionMismatch("Slater determinants of different dimensions"))

        numorbitals = length(a)

        terms = NBodyTerm[]

        ao = a.orbitals
        bo = b.orbitals

        akN = similar(ao, $N)
        blN = similar(bo, $N)

        iN = zeros(Int, $N)
        jN = zeros(Int, $N)

        aiN = similar(ao, $N)
        bjN = similar(bo, $N)

        ks,ls = nonzero_minors($N, overlap)

        for (k,l) in zip(ks,ls)
            Dkl = cofactor(k, l, overlap) # Determinantal overlap of all other orbitals
            iszero(Dkl) && continue

            for d in eachindex(k)
                akN[d] = ao[k[d]]
                blN[d] = bo[l[d]]
            end

            # Generate all distinct permutations of the N first
            # orbitals of the current total permutation.
            @above_diagonal_loop $N i N begin
                Base.Cartesian.@nexprs $N d -> iN[d] = i_d
                si = permutation_sign(iN)
                Base.Cartesian.@nexprs $N d -> aiN[d] = akN[i_d]

                @anti_diagonal_loop $N j N begin
                    Base.Cartesian.@nexprs $N d -> jN[d] = j_d
                    sj = permutation_sign(jN)
                    Base.Cartesian.@nexprs $N d -> bjN[d] = blN[j_d]

                    me = OrbitalMatrixElement(copy(aiN), op, copy(bjN))
                    iszero(me) && continue
                    nbme = convert(NBodyMatrixElement, si*sj*Dkl*me)
                    append!(terms, nbme.terms)
                end
            end
        end

        NBodyMatrixElement(terms)
    end
end

"""
    NBodyMatrixElement(a, op, b, overlap)

Generate the matrix element of `op`, a linear combination of
[`NBodyOperator`](@ref), between the Slater determinants `a` and `b`,
according to the Löwdin rules. The matrix `overlap` contains the
mutual overlaps between all single-particle orbitals in the Slater
determinants. If the orbitals are all orthogonal, the Löwdin rules
collapse to the Slater–Condon rules.
"""
function NBodyMatrixElement(a::SlaterDeterminant, op::LinearCombinationOperator, b::SlaterDeterminant, overlap)
    terms = NBodyTerm[]
    for (o,coeff) in op.operators
        iszero(coeff) && continue
        nbme = coeff*NBodyMatrixElement(a, o, b, overlap)
        !iszero(nbme) && append!(terms, nbme.terms)
    end
    NBodyMatrixElement(terms)
end

"""
    transform(f::Function, nbme::NBodyMatrixElement)

Transform integrals of the the N-body matrix element `nbme` according
to the function `f`, which should accept a single
[`NBodyTermFactor`](@ref) as its argument, and return a
[`NBodyMatrixElement`](@ref). This is useful for adapting energy
expressions to specific symmetries of the system under consideration.
"""
function transform(f::Function, nbme::NBodyMatrixElement)
    iszero(nbme) && return nbme
    sum(map(t -> transform(f, t), nbme.terms))
end

# * Overlap matrices

function overlap_matrix!(S::AbstractSparseMatrix{<:NBodyTerm}, a::SlaterDeterminant, b::SlaterDeterminant)
    length(a) == length(b) || throw(ArgumentError("Configurations not of same length: $a & $b"))
    for (i,ao) in enumerate(a.orbitals)
        for (j,bo) in enumerate(b.orbitals)
            if ao == bo
                S[i,j] = one(NBodyTerm)
            end
        end
    end

    S
end

"""
    overlap_matrix(a::SlaterDeterminant, b::SlaterDeterminant[, overlaps=[]])

Generate the single-particle orbital overlap matrix, between the
orbitals in the Slater determinants `a` and `b`. All orbitals are
assumed to be orthogonal, except for those which are given in
`overlaps`.

# Examples

First we define two Slater determinants that have some orbitals in common:
```jldoctest overlap_matrix
julia> sa = SlaterDeterminant([:i, :j, :l,:k̃])
i(1)j(2)l(3)k̃(4) - i(1)j(2)l(4)k̃(3) - i(1)j(3)l(2)k̃(4) + i(1)j(3)l(4)k̃(2) + …  + i(4)j(1)l(3)k̃(2) + i(4)j(2)l(1)k̃(3) - i(4)j(2)l(3)k̃(1) - i(4)j(3)l(1)k̃(2) + i(4)j(3)l(2)k̃(1)

julia> sb = SlaterDeterminant([:i, :j, :k, :l̃])
i(1)j(2)k(3)l̃(4) - i(1)j(2)k(4)l̃(3) - i(1)j(3)k(2)l̃(4) + i(1)j(3)k(4)l̃(2) + …  + i(4)j(1)k(3)l̃(2) + i(4)j(2)k(1)l̃(3) - i(4)j(2)k(3)l̃(1) - i(4)j(3)k(1)l̃(2) + i(4)j(3)k(2)l̃(1)
```
The orbital overlap matrix by default is
```jldoctest overlap_matrix
julia> overlap_matrix(sa, sb)
4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 2 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
```
which has only two non-zero entries, since only two of the orbitals
are common between the Slater determinants `sa` and `sb`.

We can then define that the orbitals `k̃` and `l̃` are non-orthogonal:
```jldoctest overlap_matrix
julia> overlap_matrix(sa, sb, [OrbitalOverlap(:k̃,:l̃)])
4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 3 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
  [4, 4]  =  ⟨k̃|l̃⟩
```
We can even specify that the orbital `k̃` is non-orthogonal to _itself_
(this can be useful when the `k̃` is a linear combination of orthogonal
orbitals):
```jldoctest overlap_matrix
julia> overlap_matrix(sa, sa, [OrbitalOverlap(:k̃,:k̃)])
4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 4 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
  [3, 3]  =  1
  [4, 4]  =  ⟨k̃|k̃⟩
```
Notice that this overlap matrix was calculated between the Slater determinant `sa` and itself.
"""
function overlap_matrix(a::SlaterDeterminant, b::SlaterDeterminant, overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[])
    m = length(a)
    S = spzeros(NBodyTerm, m, m)
    overlap_matrix!(S, a, b)

    for overlap in overlaps
        for (oa,ob) in [(overlap.a,overlap.b),(overlap.b,overlap.a)]
            i = findfirst(isequal(oa), a.orbitals)
            isnothing(i) && continue
            j = findfirst(isequal(ob), b.orbitals)
            isnothing(j) && continue
            S[i,j] = OrbitalOverlap(oa,ob)
        end
    end

    S
end

# * Energy expression matrix

"""
    EnergyExpression

An energy expression is given by an energy matrix, or interaction
matrix, sandwiched between a vector of mixing coefficients: `E =
c'H*c`, where `c` are the mixing coefficients and `H` the energy
matrix.
"""
const EnergyExpression = AbstractMatrix{NBodyMatrixElement}

"""
    Matrix(op::QuantumOperator, slater_determinants[, overlaps])

Generate the matrix corresponding to the quantum operator `op`,
between the different `slater_determinants`. It is possible to specify
non-orthogonalities between single-particle orbitals in `overlaps`.
"""
function Base.Matrix(op::QuantumOperator, slater_determinants::VSD,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) where {VSD<:AbstractVector{<:SlaterDeterminant}}
    m = length(slater_determinants)

    M = spzeros(NBodyMatrixElement, m, m)

    for (i,a) in enumerate(slater_determinants)
        for (j,b) in enumerate(slater_determinants)
            S = overlap_matrix(a,b,overlaps)

            me = NBodyMatrixElement(a,op,b, S)
            iszero(me) && continue
            M[i,j] = me
        end
    end

    M
end

export OrbitalOverlap, numbodies, overlap_matrix, transform, EnergyExpression, isdependent
