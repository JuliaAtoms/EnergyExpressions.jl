import AtomicLevels: AbstractOrbital

"""
    Conjugate(orbital)

Type representing the conjugation of an `orbital`.

# Examples

```jldoctest
julia> Conjugate(:a)
:a†
```
"""
struct Conjugate{O}
    orbital::O
end

function Base.show(io::IO, co::Conjugate{O}) where O
    show(io, co.orbital)
    write(io, "†")
end

"""
    conj(o::AbstractOrbital)

Convenience function to conjugate an `AbstractOrbital`.

# Examples

```jldoctest
julia> conj(o"1s")
1s†
```
"""
Base.conj(o::O) where {O<:AbstractOrbital} = Conjugate(o)

"""
    conj(o::Conjugate)

Convenience function to unconjugate a conjugated orbital.

# Examples

```jldoctest
julia> conj(Conjugate(:a))
:a
```
"""
Base.conj(co::Conjugate{O}) where O = co.orbital

export Conjugate
