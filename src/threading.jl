# https://julialang.org/blog/2023/07/PSA-dont-use-threadid/
using Base.Threads: nthreads, @spawn
function tmapreduce(f, op, itr; tasks_per_thread::Int = 2, kwargs...)
    chunk_size = max(1, length(itr) ÷ (tasks_per_thread * nthreads()))
    tasks = map(Iterators.partition(itr, chunk_size)) do chunk
        @spawn mapreduce(f, op, chunk; kwargs...)
    end
    mapreduce(fetch, op, tasks; kwargs...)
end
