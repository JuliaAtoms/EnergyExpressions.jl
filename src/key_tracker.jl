# A wrapper around a Dict that remembers at which index a certain key
# was inserted, and returns this index when queried. When asked for
# all keys, they are returned in insertion order.
struct KeyTracker{T,Data<:AbstractDict{T,Int}}
    data::Data
end

KeyTracker{T}() where T = KeyTracker(Dict{T,Int}())

Base.get!(kt::KeyTracker, i) =
    get!(kt.data, i, length(kt.data)+1)

function Base.keys(kt::KeyTracker)
    kv = collect(pairs(kt.data))
    first.(kv)[sortperm(last.(kv))]
end

struct LockedDict{K,V,Lock} <: AbstractDict{K,V}
    d::Dict{K,V}
    lock::Lock
end

LockedDict{K,V}() where {K,V} = LockedDict(Dict{K,V}(), ReentrantLock())

Base.pairs(d::LockedDict) = pairs(d.d)
Base.length(d::LockedDict) = length(d.d)

function Base.get!(d::LockedDict, i, default)
    i ∈ keys(d.d) && return d.d[i]
    lock(d.lock) do
        get!(d.d, i, default)
    end
end
