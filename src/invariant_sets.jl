"""
    coupled_states(E[; i₀=1])

Find all states coupled by the energy expression `E`, starting from
the state with index `i₀`. This can be useful to reduce the necessary
basis or to generate invariant sets for split-operator propagation.
"""
function coupled_states(E::AbstractSparseMatrix; i₀=1)
    m = size(E,1)
    visited = falses(m)
    visited[i₀] = true

    rows = rowvals(E)

    istack = [i₀]
    while !isempty(istack)
        icur = pop!(istack)
        neighbours = rows[nzrange(E, icur)]
        for n in neighbours
            if !visited[n]
                push!(istack, n)
                visited[n] = true
            end
        end
    end

    visited
end

function invariant_sets(E::AbstractSparseMatrix)
    m = size(E,1)
    visited = falses(m)

    sets = Vector{Vector{Int}}()

    while !all(visited)
        icur = findfirst(.!visited)
        set = coupled_states(E; i₀=icur)
        visited[:] .|= set
        j = findall(set)
        # If the state icur is the only one in its set, it may be that
        # it is actually not coupled to anything.
        length(j) == 1 && iszero(E[icur,j]) && continue
        push!(sets, j)
    end
    sets
end

export coupled_states, invariant_sets
