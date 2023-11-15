# * Binary utilities

function showbin(io::IO, c::BitVector)
    # Printing order = LSB ... MSB
    for (i,e) in enumerate(c)
        write(io, e ? "1" : ".")
        i%4 == 0 && write(io, " ")
    end
end
showbin(c::BitVector) = showbin(stdout, c)

showbin(io::IO, i::Integer) = showbin(io, BitVector(digits(i, base=2)))

struct BitPattern{O}
    obj::O
end

Base.show(io::IO, b::BitPattern) = showbin(io, b.obj)

macro showbin(var)
    name = string(var)
    quote
        println($name, " = ", BitPattern($(esc(var))))
    end
end

# * Orbitals

"""
    Orbitals(orbitals, overlaps, has_overlap, non_orthogonalities)

Structure storing a common set of `orbitals`, along with possible
`overlaps` between them, in case of non-orthogonalities. `has_overlap`
is a boolean matrix indicates if a pair of orbitals have overlap,
either due to non-orthogonality or if they are the same
orbital. `non_orthogonalities` is a boolean vector that indicates if a
specific orbital is non-orthogonal to _any_ other orbital in the set
of orbitals. This structure is used internally by
[`BitConfigurations`](@ref).
"""
struct Orbitals{O,Overlaps,HasOverlap,NonOrthogonalities}
    orbitals::Vector{O}
    overlaps::Overlaps
    has_overlap::HasOverlap
    non_orthogonalities::NonOrthogonalities
    function Orbitals(orbitals::Vector{O}, overlaps::AbstractVector{<:OrbitalOverlap}) where O
        orb_map = Dict(o => i for (i,o) in enumerate(orbitals))
        n = length(orbitals)
        S = spdiagm(n,n,0 => fill(one(NBodyTerm), n))
        for oo in overlaps
            i = orb_map[oo.a]
            j = orb_map[oo.b]
            S[i,j] = OrbitalOverlap(oo.a, oo.b)
            S[j,i] = OrbitalOverlap(oo.b, oo.a)
        end
        N = falses(n,n)
        for oo in overlaps
            i = orb_map[oo.a]
            j = orb_map[oo.b]
            N[i,j] = true
            N[j,i] = true
        end
        has_overlap =  N .| Matrix(I, size(N))
        non = map(|, vec(reduce(|, N, dims=1)), vec(reduce(|, N, dims=2)))
        new{O,typeof(S),BitMatrix,BitVector}(orbitals, S, has_overlap, non)
    end
    function Orbitals(orbitals::Vector{O}) where O
        no = length(orbitals)
        has_overlap = Matrix(I,no,no)
        new{O,UniformScaling{Bool},BitMatrix,Nothing}(orbitals, I, has_overlap, nothing)
    end
end

Base.length(orbitals::Orbitals) = length(orbitals.orbitals)
Base.eachindex(orbitals::Orbitals) = eachindex(orbitals.orbitals)
Base.iterate(orbitals::Orbitals, args...) = iterate(orbitals.orbitals, args...)
Base.getindex(orbitals::Orbitals, i) = orbitals.orbitals[i]
AtomicLevels.orbitals(orbitals::Orbitals) = orbitals.orbitals

Base.:(==)(a::Orbitals, b::Orbitals) =
    a.orbitals == b.orbitals &&
    a.overlaps == b.overlaps &&
    a.has_overlap == b.has_overlap &&
    a.non_orthogonalities == b.non_orthogonalities

# * BitConfiguration

struct BitConfiguration{Orbitals,Occupations}
    orbitals::Orbitals
    occupations::Occupations
end

Base.show(io::IO, c::BitConfiguration, sel=eachindex(c.orbitals)) =
    write(io, join(string.(c.orbitals[filter(∈(sel), orbital_numbers(c))]), " "))

showbin(io::IO, c::BitConfiguration) = showbin(io, c.occupations)
showbin(c::BitConfiguration) = showbin(c.occupations)

Base.:(==)(a::BitConfiguration, b::BitConfiguration) =
    a.orbitals == b.orbitals &&
    a.occupations == b.occupations

num_particles(c::BitConfiguration) = count(c.occupations)
orbital_numbers(c::BitConfiguration) = findall(c.occupations)
get_orbitals(c::BitConfiguration) = c.orbitals[orbital_numbers(c)]

function configuration_orbital_numbers(c::BitConfiguration, mask::BitVector)
    occ = copy(c.occupations)
    n = 0
    ons = Vector{Int}()
    while !iszero(occ)
        occ[1] && (n += 1)
        mask[1] && push!(ons, n)
        occ <<= 1
        mask <<= 1
    end
    ons
end

for f in (:(&), :(|), :xor)
    @eval begin
        Base.$(f)(a::BitConfiguration, b::BitConfiguration) = map($f, a.occupations, b.occupations)
        Base.$(f)(a::BitConfiguration, b::BitVector) = map($f, a.occupations, b)
        Base.$(f)(a::BitVector, b::BitConfiguration) = map($f, a, b.occupations)
    end
end

#=
Adapted from

- Scemama, A., & Giner, E. (2013). An efficient implementation of
  Slater–Condon rules. CoRR, arXiv:1311.6244.

=#

num_excitations(a::BitConfiguration, b::BitConfiguration) = count(a ⊻ b)

holes_mask(a::BitConfiguration, b::BitConfiguration) = (a ⊻ b) & a
holes(a::BitConfiguration, b::BitConfiguration) = findall(holes_mask(a, b))

particles_mask(a::BitConfiguration, b::BitConfiguration) = (a ⊻ b) & b
particles(a::BitConfiguration, b::BitConfiguration) = findall(particles_mask(a, b))

function holes_particles_mask(a::BitConfiguration, b::BitConfiguration)
    x = a ⊻ b
    x, x & a, x & b
end

function holes_particles(a::BitConfiguration, b::BitConfiguration)
    x, hm, pm = holes_particles_mask(a, b)
    findall(hm), findall(pm)
end

function relative_sign(a::BitConfiguration, b::BitConfiguration,
                       holes, particles; verbosity=0)
    isempty(holes) && return 1

    perm = 0
    av = copy(a.occupations)
    n = length(av)

    if verbosity > 0
        print("   ")
        @showbin a
        print("   ")
        @showbin b
    end

    @inbounds for (h,p) in zip(holes, particles)
        lo,hi = minmax(h,p)
        mask = vcat(falses(lo), trues(hi-lo-1), falses(n-hi+1))
        perm += count(map(&, mask, av))
        if verbosity > 0
            @show h,p,lo,hi
            ma = map(&, mask, av)
            print("  ")
            @showbin av
            @showbin mask
            print("  ")
            @showbin ma
        end
        av[h] = false
        av[p] = true
    end
    @assert av == b.occupations

    iseven(perm) ? 1 : -1
end

# * BitConfigurations

"""
    BitConfigurations(orbitals, configurations)

Represent collection of `configurations` as bit vectors, where `true`
values indicate that specific `orbitals` are occupied.

# Example

```julia-repl
julia> bcs = BitConfigurations([[:a,:b,:c], [:x,:b,:c], [:a,:y,:c], [:a,:b,:z]])
6-orbital 4-configuration BitConfigurations

1: a b c
2: a -> x
3: b -> y
4: c -> z

julia> h = FieldFreeOneBodyHamiltonian()
ĥ₀

julia> Matrix(bcs, h)
4×4 SparseMatrixCSC{NBodyMatrixElement, Int64} with 10 stored entries:
 (a|a) + (b|b) + (c|c)  (a|x)                  - (b|y)                (c|z)
 (x|a)                  (b|b) + (c|c) + (x|x)  ⋅                      ⋅
 - (y|b)                ⋅                      (a|a) + (c|c) + (y|y)  ⋅
 (z|c)                  ⋅                      ⋅                      (a|a) + (b|b) + (z|z)
```
"""
struct BitConfigurations{Orbitals}
    orbitals::Orbitals
    configurations::BitMatrix
end

get_orbitals(c::Configuration) = c.orbitals
get_orbitals(v::AbstractVector) = v
AtomicLevels.orbitals(bcs::BitConfigurations) = orbitals(bcs.orbitals)

function BitConfigurations(orbitals::Orbitals,
                           configurations::AbstractVector)
    orb_map = Dict(o => i for (i,o) in enumerate(orbitals))
    C = falses(length(orbitals), length(configurations))
    p = Progress(prod(size(C)))
    for (j,c) in enumerate(configurations)
        for o in get_orbitals(c)
            C[orb_map[o],j] = true
            ProgressMeter.next!(p)
        end
    end
    BitConfigurations(orbitals, C)
end

function BitConfigurations(configurations::AbstractVector,
                           overlaps=OrbitalOverlap[])
    orbitals = (!isempty(configurations) ?
        unique(reduce(vcat, map(get_orbitals, configurations))) :
        [])
    BitConfigurations(Orbitals(orbitals, overlaps),
                      configurations)
end

Base.getindex(m::BitConfigurations, i::Integer) = BitConfiguration(m.orbitals, m.configurations[:,i])
Base.getindex(m::BitConfigurations, i) = BitConfigurations(m.orbitals, m.configurations[:,i])
Base.length(m::BitConfigurations) = size(m.configurations, 2)
Base.isempty(m::BitConfigurations) = length(m) == 0

Base.:(==)(a::BitConfigurations, b::BitConfigurations) =
    a.orbitals == b.orbitals &&
    a.configurations == b.configurations

function AtomicLevels.core(m::BitConfigurations)
    isempty(m) && return 1:0
    # Find all initial orbitals which are occupied in all
    # configurations.
    i = something(findfirst(.!vec(reduce(&, m.configurations, dims=2))),
                  length(m.orbitals)+1) - 1
    sel = 1:i
    if i > 0 && !isnothing(m.orbitals.overlaps)
        # We require that the core orbitals are all canonical ones,
        # i.e. no non-orthogonalities are allowed. This simplifies the
        # energy expression, since the core can then only contribute
        # through Slater–Condon rules.
        non = m.orbitals.non_orthogonalities[sel]
        imax = sel[end]
        sel = 1:min(imax, something(findfirst(non), imax+1))
    end
    sel
end

function Base.show(io::IO, m::BitConfigurations)
    norb = length(m.orbitals)
    ncfg = length(m)
    write(io, "$(norb)-orbital $(ncfg)-configuration BitConfigurations")
end

function Base.show(io::IO, ::MIME"text/plain", m::BitConfigurations)
    show(io, m)
    println(io)
    c = core(m)
    !isempty(c) && write(io, "Common core: $(c)")
    isempty(m) && return
    println(io)
    ncfg = length(m)
    fmt = FormatExpr("{1:$(length(digits(ncfg)))d}: ")
    printfmt(io, fmt, 1)
    ref = m[1]
    sel = (isempty(c) || length(m) == 1 ? 0 : c[end])+1:length(m.orbitals)
    show(io, ref, sel)
    show_conf = i -> begin
        println(io)
        exc = m[i]
        h,p = holes_particles(ref, exc)
        printfmt(io, fmt, i)
        write(io, join(string.(m.orbitals[h]), " "))
        write(io, " -> ")
        write(io, join(string.(m.orbitals[p]), " "))
    end
    screen_rows = fld(max(3,first(displaysize(io))-6), 2)
    first_block = 2:min(ncfg, screen_rows)
    second_block = max(1,ncfg - screen_rows):ncfg
    if !isempty(first_block) && isempty((2:first_block[end]+2) ∩ second_block)
        foreach(show_conf, first_block)
        print(io, "\n⋮")
        foreach(show_conf, second_block)
    else
        foreach(show_conf, 2:second_block[end])
    end
end

function orbital_overlap_matrix(sd::BitConfigurations, i, j)
    a = sd[i]
    b = sd[j]

    sd.orbitals.overlaps[orbital_numbers(a), orbital_numbers(b)]
end

"""
    non_zero_cofactors(sd, N, i, j)

Find all non-zero cofactors of the orbital overlap matrix between the
Slater determinants `i` & `j` of `sd`, generated when striking out `N`
rows & columns. This routine is tailored towards the case when few
non-orthogonalities are present, e.g. approximately proportional to
the number of orbitals.

Non-orthogonality between spin-orbitals is handled by dividing the
them into two subspaces:

1. The orthogonal spin-orbitals that are common to both Slater
   determinants (core orbitals),

2. All non-orthogonal orbitals, and the orbitals which differ between
   the Slater determinants (i.e. holes of `i` and particles of `j`).

The relative phase between the Slater determinants is determined by
group `2` alone, by permuting the particles to the positions of the
holes, we find this phase. We can then formally permute them together
to a diagonal block at lower-right corner of the orbital overlap
matrix without incurring a phase change, since we need to permute the
same number of rows and columns. We thus get this structure:

                     ╷       ╷
                     │ 1 │   │
     det(Sᵢⱼ) = (-)ᵏ │───┼───│
                     │   │ 2 │
                     ╵       ╵

where `k` is decided by the permutation necessary to put the particles
in the positions of the holes.

Obviously, the determinant of the orbital matrix is now given by
`det(Sᵢⱼ) = (-)ᵏ*det(2)`, since we trivially have `det(1)==1`.

Depending on the rank of `2` (determined by the number of
hole–particle pairs and which spin-orbitals are non-orthogonal), we
need to strike out at least `size(2,1)-rank(2)` rows/columns from `2`,
and at most `min(N,size(2,1))`, i.e. for each value of `n ∈
size(2,1)-rank(2):min(N,size(2,1))`, we need to additionally strike
out `m = N - n` rows from `1`, but since the determinant of subspace
`1` is unity, regardless of how many rows/columns we've stricken out,
this is a trivial excercise. Of course, we also require that `m ≤
size(1,1)`.
"""
function non_zero_cofactors(sd::BitConfigurations, N, i, j; verbosity=0)
    verbosity > 0 && @show N
    a = sd[i]
    b = sd[j]
    np = num_particles(a)
    np == num_particles(b) ||
        throw(ArgumentError("Configurations do not have the same number of particles"))

    ks = Vector{Vector{Int}}()
    ls = Vector{Vector{Int}}()
    Ds = Vector{NBodyMatrixElement}()

    # First deal with trivial special cases.
    N > np && return ks, ls, Ds

    aon = orbital_numbers(a)
    bon = orbital_numbers(b)

    if N == np
        push!(ks, aon)
        push!(ls, bon)
        push!(Ds, one(NBodyMatrixElement))
        return ks, ls, Ds
    end

    isnothing(sd.orbitals.non_orthogonalities) &&
        error("Straight Slater–Condon not yet implemented")

    if verbosity > 0
        S = orbital_overlap_matrix(sd, i, j)
        println("Overlap matrix:")
        show(stdout, "text/plain", S)
        println()
        i == j && @show det(S)
    end

    v = sd.orbitals.non_orthogonalities
    excitations, hmask, pmask = holes_particles_mask(a, b)

    common = a & b
    # For these, we may employ Slater–Condon rules
    canonical_orbitals = findall(map(&, common, .~v))
    verbosity > 1 && @show canonical_orbitals

    # For these, we have to apply the Löwdin rules, but we can
    # precompute the determinant of this subspace once, for all those
    # determinant minors which do no strike out rows/columns in this
    # subspace.
    x = map(|, v, excitations)
    anon = findall(a & x)
    bnon = findall(b & x)

    h = findall(hmask)
    p = findall(pmask)

    S2 = sd.orbitals.overlaps[anon, bnon]
    non2 = sd.orbitals.has_overlap[anon, bnon]
    if verbosity > 0
        show(stdout, "text/plain", S2)
        println()
        show(stdout, "text/plain", sparse(non2))
        println()
        isempty(S2) && println("We're essentially in Slater–Condon land")
    end

    # Minimum order of Γ⁽ᵖ⁾
    prmin = count(iszero, reduce(|, non2, dims=2))
    pcmin = count(iszero, reduce(|, non2, dims=1))
    pmin = max(prmin, pcmin)

    if N < pmin
        verbosity > 0 && println("Too many vanishing rows/columns")
        return ks, ls, Ds
    end

    verbosity > 0 && println("We need to strike out $(pmin)..$(min(N, length(anon))) rows/columns from S2")

    rs = relative_sign(a, b, h, p, verbosity=verbosity)
    verbosity > 0 && @show rs

    if N == 0
        # As a special case, for the zero-body operator, return the
        # properly phased determinant of S2.
        push!(ks, [])
        push!(ls, [])
        push!(Ds, rs*det(S2))
        return ks, ls, Ds
    end

    for n = pmin:min(N, length(anon))
        # n is how many rows/columns we're striking out from S2;
        # consequently m is how many rows/columns we have to strike
        # out from S1
        m = N - n
        verbosity > 0 && @show n, m
        kk,ll = nonzero_minors(n, S2)
        for (k,l) in zip(kk,ll)
            verbosity > 0 && @show k,l
            Dkl = cofactor(k, l, S2)
            iszero(Dkl) && continue
            verbosity > 0 && @show Dkl
            # We now employ the Slater–Condon rules for S1
            for co in combinations(canonical_orbitals, m)
                ko = vcat(anon[k], co)
                lo = vcat(bnon[l], co)
                push!(ks, ko)
                push!(ls, lo)
                push!(Ds, rs*Dkl)
            end
        end
    end

    ks, ls, Ds
end

function NBodyMatrixElement(bcs::BitConfigurations, op::NBodyOperator{N}, i, j) where N
    orbitals = bcs.orbitals.orbitals
    ks,ls,Ds = non_zero_cofactors(bcs, N, i, j)
    NBodyMatrixElement(orbitals, op, orbitals,
                       zip(ks, ls, Ds))
end

function NBodyMatrixElement(bcs::BitConfigurations, op::LinearCombinationOperator, i, j)
    terms = NBodyTerm[]
    for (o,coeff) in op.operators
        iszero(coeff) && continue
        nbme = coeff*NBodyMatrixElement(bcs, o, i, j)
        !iszero(nbme) && append!(terms, nbme.terms)
    end
    NBodyMatrixElement(terms)
end

function Base.Matrix(bcs::BitConfigurations, rows, op::QuantumOperator, cols; verbosity=0, kwargs...)
    m,n = length(rows), length(cols)

    # https://discourse.julialang.org/t/threading-race-condition-when-pushing-to-arrays/99567/3
    # Allocate channels for m rows.
    Is = Channel{Vector{Int}}(m)
    Js = Channel{Vector{Int}}(m)
    Vs = Channel{Vector{NBodyMatrixElement}}(m)

    p = if verbosity > 0
        @info "Generating energy expression" op
        Progress(m*n)
    end

    Threads.@threads for i in eachindex(rows)
        r = rows[i]
        I = Int[]
        J = Int[]
        V = NBodyMatrixElement[]
        for (j,c) in enumerate(cols)
            me = NBodyMatrixElement(bcs, op, r, c)

            !isnothing(p) && ProgressMeter.next!(p)

            iszero(me) && continue

            push!(I, i)
            push!(J, j)
            push!(V, me)
        end
        put!(Is, I)
        put!(Js, J)
        put!(Vs, V)
    end
    close(Is); close(Js); close(Vs)

    sparse(reduce(vcat, Is),
           reduce(vcat, Js),
           reduce(vcat, Vs), m, n)
end

function Base.Matrix(bcs::BitConfigurations, op::QuantumOperator; left=Colon(), right=Colon(), kwargs...)
    ms = 1:length(bcs)
    Matrix(bcs, ms[left], op, ms[right]; kwargs...)
end

export BitConfigurations
