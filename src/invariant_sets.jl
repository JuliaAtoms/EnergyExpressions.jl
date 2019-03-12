"""
    coupled_states(E[; i₀=1])

Find all states coupled by the energy expression `E`, starting from
the state with index `i₀`. This can be useful to reduce the necessary
basis or to generate invariant sets for split-operator propagation.
"""
function coupled_states(E::EM; i₀=1) where {EM<:EnergyExpression}
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

function invariant_sets(H::QuantumOperator)
    m = size(hamiltonians[1],1)
    visited = falses(m)

    sets = Vector{Vector{Bool}}()

    while !all(visited)
        icur = findfirst(.!visited)
        set = coupled_states(H; i₀=icur)
        push!(sets, set)
        visited[:] .|= set
    end
    sets
end

export coupled_states, invariant_sets
