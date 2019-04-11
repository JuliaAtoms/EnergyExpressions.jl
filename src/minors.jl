count_nonzeros(IJ, mn) = [count(isequal(ij), IJ) for ij in 1:mn]

"""
    nonzero_minors(N, overlap) -> (ks,ls)

Find all minor determinants of order `N` of the orbital `overlap`
matrix that do not vanish, i.e. all non-vanishing minors are
guaranteed to be present, but not all of the returned minors are
guaranteed to be non-zero. Vanishing minors returned arise when the
overlap matrix is rank deficient, which is unlikely to happen when
computing energy expressions, but must still be guarded against. This
is most easily checked by actually calculating the [`cofactor`](@ref),
which is most likely desired anyway.
"""
function nonzero_minors(N::Integer, overlap::AbstractSparseMatrix)
    ks = Vector{Vector{Int}}()
    ls = Vector{Vector{Int}}()

    (I,J,V) = findnz(overlap)
    m,n = size(overlap)
    m == n || throw(ArgumentError("Overlap matrix must be square"))

    # First deal with trivial special cases.
    N > m && return ks, ls
    if N == m
        push!(ks, collect(1:m))
        push!(ls, collect(1:m))
        return ks, ls
    end

    # Number of non-zero matrix elements for each row/column.
    nzrowcount = count_nonzeros(I, m)
    nzcolcount = count_nonzeros(J, n)
    # Vanishing rows/columns.
    zrows = findall(iszero, nzrowcount)
    zcols = findall(iszero, nzcolcount)
    # Minimum order of Γ⁽ᵖ⁾
    prmin = length(zrows)
    pcmin = length(zcols)
    pmin = max(prmin, pcmin)

    # Then deal with next set of trivial special cases.
    N < pmin && return ks, ls
    if N == prmin && N == pcmin
        push!(ks, zrows)
        push!(ls, zcols)
        return ks, ls
    end

    # The general case.
    Nr = N - prmin
    Nc = N - pcmin
    # Find all rows supporting a specific column j.
    I′ = [I[findall(isequal(j), J)] for j in 1:n]
    # Find all columns supporting a specific row i.
    J′ = [J[findall(isequal(i), I)] for i in 1:m]

    for k in combinations(unique(I),Nr)
        # Find all columns supported by any of the rows in the
        # combination k.
        cols = unique(vcat([J′[i] for i in k]...))
        # For each combination k of rows, we can at most strike
        # out Nc columns.
        colbudget = Nc
        colsmuststrike = Vector{Int}()

        for j in cols
            # Test if column j is solely supported by the row
            # combination k.
            colsupport = setdiff(I′[j], k)
            if isempty(colsupport)
                colbudget -= 1
                colbudget < 0 && break
                push!(colsmuststrike, j)
            end
        end
        colbudget < 0 && continue

        Ncm = length(colsmuststrike)

        # The columns we must strike out fill the column budget.
        if Ncm == Nc
            push!(ks, sort(vcat(k,zrows)))
            push!(ls, sort(vcat(colsmuststrike,zcols)))
            continue
        end

        # Find all combinations of the remaining columns.
        colchoices = combinations(setdiff(unique(J), colsmuststrike), Nc-Ncm)
        for l in colchoices
            # Find all rows that would be affected if column
            # combination l is stricken out.
            colsupport = unique(vcat([I′[j] for j in l]...))
            supported = true
            for i in colsupport
                # If the row i supporting one of the columns in the
                # column combination l is already a candidate to be
                # stricken out, it does not matter if its support
                # vanishes.
                i ∈ k && continue
                # Otherwise, if its support vanishes, then l is an
                # unviable column combination.
                rowsupport = setdiff(J′[i], l)
                if isempty(rowsupport)
                    supported = false
                    break
                end
            end
            if supported
                push!(ks, sort(vcat(k,zrows)))
                push!(ls, sort(vcat(l,colsmuststrike,zcols)))
            end
        end
    end

    ks, ls
end
