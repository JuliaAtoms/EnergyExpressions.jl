function compare_vectors(a, b)
    na = length(a)
    nb = length(b)
    # We don't bother with the case when one element in a is formally the
    # same as the sum of multiple elements in b.
    na == nb == 0 && return true

    # We need to keep track of which elements with compared to, such
    # that we use the same element in b multiple time to match
    # different elements in a.
    is = Vector(1:nb)
    test_element(t) = findfirst(i -> t == b[i], is)
    for e in a
        iszero(e) && continue
        i = test_element(e)
        isnothing(i) && return false
        deleteat!(is, i)
    end
    isempty(is) || all(i -> iszero(b[i]), is)
end
