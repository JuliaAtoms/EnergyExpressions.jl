"""
    above_diagonal_loop(N, itersym, imax, args...)

Generate `N` Cartesian loops for the iteration variables
`itersym_{1:N}`, where `itersym_N ∈ 1:imax`, `itersym_{N-1} ∈
itersym_N+1:imax`, etc, i.e. above the hyper-diagonal of the
`N`-dimensional hypercube with the side `imax`. `args...` is passed on
to `Base.Cartesian._nloops`. `above_diagonal_loop` is nestable.

# Examples

```jldoctest
julia> @above_diagonal_loop 2 i 3 begin
           println("==================================")
           println("i = ", Base.Cartesian.@ntuple 2 i)
           @above_diagonal_loop 2 j 3 begin
               println("j = ", Base.Cartesian.@ntuple 2 j)
           end
       end
==================================
i = (2, 1)
j = (2, 1)
j = (3, 1)
j = (3, 2)
==================================
i = (3, 1)
j = (2, 1)
j = (3, 1)
j = (3, 2)
==================================
i = (3, 2)
j = (2, 1)
j = (3, 1)
j = (3, 2)
```
"""
macro above_diagonal_loop(N, itersym, imax, args...)
    lim = Expr(:call, :(:), Expr(:call, :(+), Expr(:curly, Symbol("$(itersym)_"), Expr(:call, :(+), :d, 1)), 1), imax)
    rng = Expr(:(->), :d, Expr(:block, Expr(:if, :(d == $N), :(1:$imax), lim)))
    Base.Cartesian._nloops(N, itersym, rng, args...)
end

"""
    anti_diagonal_loop(N, itersym, imax, args...)

Generate `N` Cartesian loops for the iteration variables
`itersym_{1:N}`, where `itersym_N ∈ 1:imax`, `itersym_{N-1} ∈
1:imax\\itersym_N`, etc, i.e. no two iteration variables have the same
values simultaneously. `args...` is passed on to
`Base.Cartesian._nloops`; however, `preexpr` is already used to skip
the diagonal elements. `anti_diagonal_loop` is nestable.

# Examples

```jldoctest
julia> @anti_diagonal_loop 3 i 3 begin
           println("-----------------------------")
           t = (Base.Cartesian.@ntuple 3 i)
           println("\$t: ", allunique(t))
           @anti_diagonal_loop 2 j 2 begin
               u = (Base.Cartesian.@ntuple 2 j)
               println("\$u: ", allunique(u))
           end
       end
-----------------------------
(3, 2, 1): true
(2, 1): true
(1, 2): true
-----------------------------
(2, 3, 1): true
(2, 1): true
(1, 2): true
-----------------------------
(3, 1, 2): true
(2, 1): true
(1, 2): true
-----------------------------
(1, 3, 2): true
(2, 1): true
(1, 2): true
-----------------------------
(2, 1, 3): true
(2, 1): true
(1, 2): true
-----------------------------
(1, 2, 3): true
(2, 1): true
(1, 2): true
```
"""
macro anti_diagonal_loop(N, itersym, imax, args...)
    rng = :(d -> 1:$imax)
    # The preexpr body generates tests for inner loops that are true
    # if the inner loop variable equals any of the outer ones; then
    # that iteration is skipped.
    prebody = Expr(:(->), :d,
                   Expr(:call, :(==),
                        Symbol("$(itersym)_e"),
                        Expr(:curly, Symbol("$(itersym)_"), :($N - d + 1))))
    pre = :(e -> (Base.Cartesian.@nany $N-e $prebody) && continue)

    Base.Cartesian._nloops(N, itersym, rng, pre, args...)
end
