# https://julialang.org/blog/2023/07/PSA-dont-use-threadid/
using Base.Threads: nthreads, @spawn
function tmapreduce(f, op, itr; tasks_per_thread::Int = 2, kwargs...)
    chunk_size = max(1, length(itr) ÷ (tasks_per_thread * nthreads()))
    tasks = map(Iterators.partition(itr, chunk_size)) do chunk
        @spawn mapreduce(f, op, chunk; kwargs...)
    end
    mapreduce(fetch, op, tasks; kwargs...)
end

function threaded_sparse_builder(fun::Function, ::Type{T}, M, N;
                                 tasks_per_thread::Int = 2,
                                 task_local::Function=()->(),
                                 verbosity=1) where T
    ci = CartesianIndices((M,N))
    nt = Threads.nthreads()
    chunk_size = max(1, length(ci) ÷ (tasks_per_thread * nt))
    chunks = Iterators.partition(ci, chunk_size)
    nc = length(chunks)

    verbosity > 1 && @info "Threaded $(M)×$(N) SparseMatrixCSC{$T} builder; $(nt) threads, $(tasks_per_thread) tasks/thread ⟹ $(nc) chunks of size $(chunk_size)"

    # https://discourse.julialang.org/t/threading-race-condition-when-pushing-to-arrays/99567/3
    # Allocate channels for nc chunks
    Is = Channel{Vector{Int}}(nc)
    Js = Channel{Vector{Int}}(nc)
    Vs = Channel{Vector{T}}(nc)

    p = verbosity > 0 ? Progress(M*N) : nothing

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            I = Int[]
            J = Int[]
            V = T[]
            tloc = task_local()

            for ii in chunk
                i,j = ii[1],ii[2]
                v = fun(i,j,tloc)

                isnothing(p) || ProgressMeter.next!(p)
                iszero(v) && continue

                push!(I, i)
                push!(J, j)
                push!(V, v)
            end
            put!(Is, I)
            put!(Js, J)
            put!(Vs, V)
        end
    end
    wait.(tasks)
    isnothing(p) || ProgressMeter.finish!(p)
    close(Is); close(Js); close(Vs)

    sparse(reduce(vcat, Is),
           reduce(vcat, Js),
           reduce(vcat, Vs), M, N)
end



function threaded_sparse_map(fun::Function, ::Type{T}, A::SparseMatrixCSC;
                             tasks_per_thread::Int = 2,
                             task_local::Function=()->(),
                             verbosity=1) where T
    M,N = size(A)
    _nnz = nnz(A)
    nt = Threads.nthreads()
    # We only parallelize over columns for now, since findnz copies
    # data, whereas nzrange just returns a range of which structural
    # non-zeros belong to a given column.
    chunk_size = max(1, N ÷ (tasks_per_thread * nt))
    chunks = Iterators.partition(1:N, chunk_size)
    nc = length(chunks)

    verbosity > 1 && @info "Threaded $(M)×$(N) SparseMatrixCSC{$T} builder, max $(_nnz) structural non-zeros; $(nt) threads, $(tasks_per_thread) tasks/thread ⟹ $(nc) chunks of size $(chunk_size)"

    # https://discourse.julialang.org/t/threading-race-condition-when-pushing-to-arrays/99567/3
    # Allocate channels for nc chunks
    Is = Channel{Vector{Int}}(nc)
    Js = Channel{Vector{Int}}(nc)
    Vs = Channel{Vector{T}}(nc)

    p = verbosity > 0 ? Progress(_nnz) : nothing

    rows = rowvals(A)
    vals = nonzeros(A)

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            I = Int[]
            J = Int[]
            V = T[]
            tloc = task_local()

            for j in chunk
                for ind in nzrange(A, j)
                    i = rows[ind]
                    v = fun(vals[ind], tloc)

                    isnothing(p) || ProgressMeter.next!(p)
                    iszero(v) && continue

                    push!(I, i)
                    push!(J, j)
                    push!(V, v)
                end
            end
            put!(Is, I)
            put!(Js, J)
            put!(Vs, V)
        end
    end
    wait.(tasks)
    isnothing(p) || ProgressMeter.finish!(p)
    close(Is); close(Js); close(Vs)

    sparse(reduce(vcat, Is),
           reduce(vcat, Js),
           reduce(vcat, Vs), M, N)
end
