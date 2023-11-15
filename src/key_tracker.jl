struct LockedDict{K,V,Lock} <: AbstractDict{K,V}
    d::Dict{K,V}
    lock::Lock
end

LockedDict{K,V}() where {K,V} = LockedDict(Dict{K,V}(), ReentrantLock())

Base.pairs(d::LockedDict) =
    lock(d.lock) do
        pairs(d.d)
    end
Base.length(d::LockedDict) =
    lock(d.lock) do
        length(d.d)
    end

Base.isempty(d::LockedDict) =
    lock(d.lock) do
        isempty(d.d)
    end
Base.iterate(d::LockedDict, args...) =
    lock(d.lock) do
        iterate(d.d, args...)
    end

function Base.get!(d::LockedDict, i, default::Function)
    lock(d.lock) do
        i ∈ keys(d.d) && return d.d[i]
        get!(d.d, i, default())
    end
end


# A wrapper around a Dict that remembers at which index a certain key
# was inserted, and returns this index when queried. When asked for
# all keys, they are returned in insertion order. Should probably be
# deprecated in favour of OrderedDict at some point.
struct KeyTracker{T,Data<:AbstractDict{T,Int}}
    data::Data
end

KeyTracker{T}() where T = KeyTracker(Dict{T,Int}())
function KeyTracker(elements::AbstractVector{T}) where T
    kt = KeyTracker{T}()
    for e in elements
        get!(kt, e)
    end
    kt
end

Base.get!(kt::KeyTracker, i) =
    get!(kt.data, i, length(kt.data)+1)

Base.get!(kt::KeyTracker{<:Any,<:LockedDict}, i) =
    get!(kt.data, i, () -> length(kt.data)+1)

Base.getindex(kt::KeyTracker, i) =
    getindex(kt.data, i)

function Base.keys(kt::KeyTracker)
    kv = collect(pairs(kt.data))
    first.(kv)[sortperm(last.(kv))]
end
