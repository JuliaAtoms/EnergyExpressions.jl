@testset "Show binary values" begin
    @test string(BitPattern(1)) == "1"
    @test string(BitPattern(2)) == ".1"
    @test string(BitPattern(3)) == "11"
    @test string(BitPattern(5)) == "1.1"
    @test string(BitPattern(16)) == ".... 1"

    @test "...1 1111 1" == @capture_out begin
        bv = BitVector(digits(504, base=2))
        showbin(bv)
    end
end

@testset "Bit configurations" begin
    @testset "Orbitals" begin
        @testset "Orthogonal orbitals" begin
            o = Orbitals(collect(1:10))
            @test o.overlaps == I

            @test size(o.has_overlap) == (10,10)
            @test o.has_overlap == I

            @test isnothing(o.non_orthogonalities)

            @test length(o) == 10
            @test eachindex(o) == 1:10
            @test o[4] == 4
        end
        @testset "Non-orthogonal orbitals" begin
            o = Orbitals(collect(1:10), [ov(3,8)])
            @test o.overlaps isa SparseMatrixCSC
            @test size(o.overlaps) == (10,10)
            @test o.overlaps[3,8].factors[1] == ov(3,8)
            @test o.overlaps[8,3].factors[1] == ov(8,3)

            @test size(o.has_overlap) == (10,10)
            @test o.has_overlap[3,8]
            @test o.has_overlap[8,3]
            @test sum(o.has_overlap - I) == 2

            @test o.non_orthogonalities == [0,0,1,0,0,0,0,1,0,0]

            @test length(o) == 10
            @test eachindex(o) == 1:10
            @test o[4] == 4
        end
    end

    @testset "BitConfiguration" begin
        orbitals = ["a", "b", "c", "d", "e", "f"]
        a = BitConfiguration(orbitals, BitVector([1,0,1]))
        b = BitConfiguration(orbitals, BitVector([1,1,0,0,1]))
        c = BitConfiguration(orbitals, BitVector([1,1,0,0,1,0]))
        d = BitConfiguration(orbitals, BitVector([1,0,1,0,1]))
        e = BitConfiguration(orbitals, BitVector([0,1,1,0,1]))
        @test b == BitConfiguration(orbitals, BitVector([1,1,0,0,1]))
        @test b ≠ c

        @test string(a) == "a c"
        @test string(b) == "a b e"

        let io = IOBuffer()
            showbin(io, a)
            @test String(take!(io)) == "1.1"
            showbin(io, b)
            @test String(take!(io)) == "11.. 1"
        end
        @test "1.1" == @capture_out begin
            showbin(a)
        end

        @test num_particles(a) == 2
        @test num_particles(b) == 3

        @test orbital_numbers(a) == [1,3]
        @test orbital_numbers(b) == [1,2,5]

        @test get_orbitals(a) == ["a", "c"]
        @test get_orbitals(b) == ["a", "b", "e"]
        @test get_orbitals(["a", "c"]) == ["a", "c"]
        @test get_orbitals(c"1s 2p") == [o"1s", o"2p"]

        @test b & d == BitVector([1,0,0,0,1])
        @test b & d.occupations == BitVector([1,0,0,0,1])
        @test b.occupations & d == BitVector([1,0,0,0,1])

        @test b | d == BitVector([1,1,1,0,1])
        @test b | d.occupations == BitVector([1,1,1,0,1])
        @test b.occupations | d == BitVector([1,1,1,0,1])

        @test b ⊻ d == BitVector([0,1,1,0,0])
        @test b ⊻ d.occupations == BitVector([0,1,1,0,0])
        @test b.occupations ⊻ d == BitVector([0,1,1,0,0])

        @test num_excitations(b,b) == 0
        @test num_excitations(b,d) == 2

        @test holes_mask(b, d) == BitVector([0,1,0,0,0])
        @test holes(b, d) == [2]

        @test particles_mask(b, d) == BitVector([0,0,1,0,0])
        @test particles(b, d) == [3]

        @test holes_particles_mask(b, d) == (BitVector([0,1,1,0,0]),BitVector([0,1,0,0,0]),BitVector([0,0,1,0,0]))
        @test holes_particles(b, d) == ([2], [3])

        test_value_and_output(1,
                              """
   a = 11.. 1
   b = 1.1. 1
(h, p, lo, hi) = (2, 3, 2, 3)
  av = 11.. 1
mask = .... .
  ma = .... .
""") do
    relative_sign(b,d,holes_particles(b,d)...,verbosity=Inf)
end
        test_value_and_output(-1,
                              """
   a = 11.. 1
   b = .11. 1
(h, p, lo, hi) = (1, 3, 1, 3)
  av = 11.. 1
mask = .1.. .
  ma = .1.. .
""") do
    relative_sign(b,e,holes_particles(b,e)...,verbosity=Inf)
end
    end

    @testset "BitConfigurations" begin
        @testset "Empty" begin
            bcs = BitConfigurations([])
            @test length(bcs) == 0
            @test isempty(bcs)
        end
        @testset "Orthogonal orbitals" begin
            bcs = BitConfigurations([[:a, :b, :c], [:a, :d, :e]])
            @test orbitals(bcs) == [:a, :b, :c, :d, :e]

            @test length(bcs) == 2
            @test !isempty(bcs)

            os = bcs.orbitals
            @test length(os) == 5
            @test os.orbitals == [:a, :b, :c, :d, :e]
            @test os.overlaps == I
            @test size(os.overlaps) == (5,5)
            @test os.has_overlap == I
            @test size(os.has_overlap) == (5,5)
            @test all(iszero, os.non_orthogonalities)

            @test bcs[1].occupations == BitVector([1,1,1,0,0])
            @test bcs[2].occupations == BitVector([1,0,0,1,1])

            @test bcs[1:2] == bcs

            @test core(bcs) == 1:1

            @test string(bcs) == "5-orbital 2-configuration BitConfigurations"

            o = one(NBodyTerm)
            z = zero(NBodyTerm)
            @test orbital_overlap_matrix(bcs, 1, 1) == [o z z
                                                        z o z
                                                        z z o]
            @test orbital_overlap_matrix(bcs, 1, 2) == [o z z
                                                        z z z
                                                        z z z]
        end

        @testset "Non-orthogonal orbitals" begin
            bd = ov(:b, :d)
            bcs = BitConfigurations([[:a, :b, :c], [:a, :d, :e]], [bd])

            @test length(bcs) == 2
            @test !isempty(bcs)

            o = one(NBodyTerm)
            z = zero(NBodyTerm)

            os = bcs.orbitals
            @test length(os) == 5
            @test os.orbitals == [:a, :b, :c, :d, :e]
            @test os.overlaps == NBodyTerm[o   z z  z z
                                                             z   o z bd z
                                                             z   z o  z z
                                                             z bd' z  o z
                                                             z   z z  z o]
            @test os.has_overlap == [1 0 0 0 0
                                     0 1 0 1 0
                                     0 0 1 0 0
                                     0 1 0 1 0
                                     0 0 0 0 1]
            @test os.non_orthogonalities == BitVector([0,1,0,1,0])

            @test bcs[1].occupations == BitVector([1,1,1,0,0])
            @test bcs[2].occupations == BitVector([1,0,0,1,1])

            @test bcs[1:2] == bcs

            @test core(bcs) == 1:1

            @test string(bcs) == "5-orbital 2-configuration BitConfigurations"

            @test orbital_overlap_matrix(bcs, 1, 1) == NBodyTerm[o z z
                                                                                   z o z
                                                                                   z z o]
            @test orbital_overlap_matrix(bcs, 1, 2) == NBodyTerm[o  z z
                                                                                   z bd z
                                                                                   z  z z]
        end

        @testset "Pretty-printing" begin
            @test strdisplay(BitConfigurations([[:a, :b, :c], [:a, :d, :e]])) == """
5-orbital 2-configuration BitConfigurations
Common core: 1:1
1: b c
2: b c -> d e"""

            @test strdisplay(BitConfigurations([[:a, :b, :c]])) == """
3-orbital 1-configuration BitConfigurations
Common core: 1:3
1: a b c"""

            @test strdisplay(BitConfigurations([[:a, :b, :c], [:d, :e]])) == """
5-orbital 2-configuration BitConfigurations

1: a b c
2: a b c -> d e"""

            @test strdisplay(BitConfigurations([])) == """
0-orbital 0-configuration BitConfigurations
"""

            @test strdisplay(BitConfigurations([[0,i] for i = 1:100])) == """
101-orbital 100-configuration BitConfigurations
Common core: 1:1
  1: 1
  2: 1 -> 2
  3: 1 -> 3
  4: 1 -> 4
  5: 1 -> 5
  6: 1 -> 6
  7: 1 -> 7
  8: 1 -> 8
  9: 1 -> 9
⋮
 91: 1 -> 91
 92: 1 -> 92
 93: 1 -> 93
 94: 1 -> 94
 95: 1 -> 95
 96: 1 -> 96
 97: 1 -> 97
 98: 1 -> 98
 99: 1 -> 99
100: 1 -> 100"""
        end
    end

    @testset "Non-zero cofactors" begin
        o = one(NBodyMatrixElement)
        z = zero(NBodyMatrixElement)

        @testset "N-body interaction, N > n" begin
            bcs = BitConfigurations([[1],[2]])
            ks, ls, Ds = non_zero_cofactors(bcs, 2, 1, 2)
            @test isempty(ks)
            @test isempty(ls)
            @test isempty(Ds)
        end

        @testset "N-body interaction, N == n" begin
            bcs = BitConfigurations([[1],[2]])
            ks, ls, Ds = non_zero_cofactors(bcs, 1, 1, 2)
            @test ks == [[1]]
            @test ls == [[2]]
            @test Ds == [o]
        end

        @testset "Pure Slater–Condon not yet implemented" begin
            os = Orbitals([1,2,3,4])
            bcs = BitConfigurations(os, [[1,2],[3,4]])
            @test_throws ErrorException non_zero_cofactors(bcs, 1, 1, 2)
        end

        @testset "Orthogonal orbitals" begin
            bcs = BitConfigurations([[1,2,3],[1,3,4],[1,2,5],[1,5,6]])

            @testset "Zero-body operators" begin
                test_value_and_output(([[]], [[]], [o]),
                                      """
N = 0
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 3 stored entries:
 1     ⋅    ⋅
  ⋅   1     ⋅
  ⋅    ⋅   1
det(S) = 1
canonical_orbitals = [1, 2, 3]
0×0 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries
0×0 SparseMatrixCSC{Bool, Int64} with 0 stored entries
We're essentially in Slater–Condon land
We need to strike out 0..0 rows/columns from S2
rs = 1
""") do
    non_zero_cofactors(bcs, 0, 1, 1, verbosity=Inf)
end
            end

            @testset "One-body operators" begin
                test_value_and_output(([[1], [2], [3]], [[1], [2], [3]], [o, o, o]),
                                      """
N = 1
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 3 stored entries:
 1     ⋅    ⋅
  ⋅   1     ⋅
  ⋅    ⋅   1
det(S) = 1
canonical_orbitals = [1, 2, 3]
0×0 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries
0×0 SparseMatrixCSC{Bool, Int64} with 0 stored entries
We're essentially in Slater–Condon land
We need to strike out 0..0 rows/columns from S2
rs = 1
(n, m) = (0, 1)
(k, l) = (Int64[], Int64[])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 1, 1, 1, verbosity=Inf)
end

                test_value_and_output(([[2]], [[4]], [-o]),
                                      """
N = 1
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 2 stored entries:
 1     ⋅    ⋅
  ⋅    ⋅    ⋅
  ⋅   1     ⋅
canonical_orbitals = [1, 3]
1×1 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries:
  ⋅
1×1 SparseMatrixCSC{Bool, Int64} with 0 stored entries:
 ⋅
We need to strike out 1..1 rows/columns from S2
   a = 111. ..
   b = 1.11 ..
(h, p, lo, hi) = (2, 4, 2, 4)
  av = 111. ..
mask = ..1. ..
  ma = ..1. ..
rs = -1
(n, m) = (1, 0)
(k, l) = ([1], [1])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 1, 1, 2, verbosity=Inf)

end

                test_value_and_output(([], [], []),
                                      """
N = 1
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 1 stored entry:
 1     ⋅    ⋅
  ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅
canonical_orbitals = [1]
2×2 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries:
  ⋅    ⋅
  ⋅    ⋅
2×2 SparseMatrixCSC{Bool, Int64} with 0 stored entries:
 ⋅  ⋅
 ⋅  ⋅
Too many vanishing rows/columns
""") do
    non_zero_cofactors(bcs, 1, 1, 4, verbosity=Inf)
end
            end

            @testset "Two-body operators" begin
                test_value_and_output(([[2, 3]], [[5, 6]], [o]),
                                      """
N = 2
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 1 stored entry:
 1     ⋅    ⋅
  ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅
canonical_orbitals = [1]
2×2 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries:
  ⋅    ⋅
  ⋅    ⋅
2×2 SparseMatrixCSC{Bool, Int64} with 0 stored entries:
 ⋅  ⋅
 ⋅  ⋅
We need to strike out 2..2 rows/columns from S2
   a = 111. ..
   b = 1... 11
(h, p, lo, hi) = (2, 5, 2, 5)
  av = 111. ..
mask = ..11 ..
  ma = ..1. ..
(h, p, lo, hi) = (3, 6, 3, 6)
  av = 1.1. 1.
mask = ...1 1.
  ma = .... 1.
rs = 1
(n, m) = (2, 0)
(k, l) = ([1, 2], [1, 2])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 2, 1, 4, verbosity=Inf)
end

                test_value_and_output(([[1, 2], [1, 3], [2, 3]], [[1, 2], [1, 3], [2, 3]], [o, o, o]),
                                      """
N = 2
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 3 stored entries:
 1     ⋅    ⋅
  ⋅   1     ⋅
  ⋅    ⋅   1
det(S) = 1
canonical_orbitals = [1, 2, 3]
0×0 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries
0×0 SparseMatrixCSC{Bool, Int64} with 0 stored entries
We're essentially in Slater–Condon land
We need to strike out 0..0 rows/columns from S2
rs = 1
(n, m) = (0, 2)
(k, l) = (Int64[], Int64[])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 2, 1, 1, verbosity=Inf)

end

                test_value_and_output(([[2, 1], [2, 3]], [[4, 1], [4, 3]], [-o, -o]),
                                      """
N = 2
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 2 stored entries:
 1     ⋅    ⋅
  ⋅    ⋅    ⋅
  ⋅   1     ⋅
canonical_orbitals = [1, 3]
1×1 SparseMatrixCSC{NBodyTerm, Int64} with 0 stored entries:
  ⋅
1×1 SparseMatrixCSC{Bool, Int64} with 0 stored entries:
 ⋅
We need to strike out 1..1 rows/columns from S2
   a = 111. ..
   b = 1.11 ..
(h, p, lo, hi) = (2, 4, 2, 4)
  av = 111. ..
mask = ..1. ..
  ma = ..1. ..
rs = -1
(n, m) = (1, 1)
(k, l) = ([1], [1])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 2, 1, 2, verbosity=Inf)
end
            end
        end

        @testset "Non-orthogonal orbitals" begin
            s = ov(2,5)
            bcs = BitConfigurations([[1,2,3],[1,3,4],[1,2,5],[1,5,6]], [s])

            test_value_and_output(([[3]], [[6]], NBodyMatrixElement[1s]),
                                  """
N = 1
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 2 stored entries:
 1     ⋅      ⋅
  ⋅   ⟨2|5⟩   ⋅
  ⋅    ⋅      ⋅
canonical_orbitals = [1]
2×2 SparseMatrixCSC{NBodyTerm, Int64} with 1 stored entry:
 ⟨2|5⟩   ⋅
  ⋅      ⋅
2×2 SparseMatrixCSC{Bool, Int64} with 1 stored entry:
 1  ⋅
 ⋅  ⋅
We need to strike out 1..1 rows/columns from S2
   a = 111. ..
   b = 1... 11
(h, p, lo, hi) = (2, 5, 2, 5)
  av = 111. ..
mask = ..11 ..
  ma = ..1. ..
(h, p, lo, hi) = (3, 6, 3, 6)
  av = 1.1. 1.
mask = ...1 1.
  ma = .... 1.
rs = 1
(n, m) = (1, 0)
(k, l) = ([2], [2])
Dkl = ⟨2|5⟩
""") do
    non_zero_cofactors(bcs, 1, 1, 4, verbosity=Inf)
end

            test_value_and_output(([[3, 1], [2, 3]], [[6, 1], [5, 6]], NBodyMatrixElement[1s, o]),
                                  """
N = 2
Overlap matrix:
3×3 SparseMatrixCSC{NBodyTerm, Int64} with 2 stored entries:
 1     ⋅      ⋅
  ⋅   ⟨2|5⟩   ⋅
  ⋅    ⋅      ⋅
canonical_orbitals = [1]
2×2 SparseMatrixCSC{NBodyTerm, Int64} with 1 stored entry:
 ⟨2|5⟩   ⋅
  ⋅      ⋅
2×2 SparseMatrixCSC{Bool, Int64} with 1 stored entry:
 1  ⋅
 ⋅  ⋅
We need to strike out 1..2 rows/columns from S2
   a = 111. ..
   b = 1... 11
(h, p, lo, hi) = (2, 5, 2, 5)
  av = 111. ..
mask = ..11 ..
  ma = ..1. ..
(h, p, lo, hi) = (3, 6, 3, 6)
  av = 1.1. 1.
mask = ...1 1.
  ma = .... 1.
rs = 1
(n, m) = (1, 1)
(k, l) = ([2], [2])
Dkl = ⟨2|5⟩
(n, m) = (2, 0)
(k, l) = ([1, 2], [1, 2])
Dkl = 1
""") do
    non_zero_cofactors(bcs, 2, 1, 4, verbosity=Inf)
end
        end
    end

    @testset "NBodyMatrixElement" begin
        s = ov(2,5)
        bcs = BitConfigurations([[1,2,3],[1,3,4],[1,2,5],[1,5,6],[1,3,5]], [s])

        I₀ = IdentityOperator{0}()
        h = FieldFreeOneBodyHamiltonian()
        g = CoulombInteraction()

        H = h + g

        @test NBodyMatrixElement(bcs, I₀, 1, 5) == -ome([], I₀, [])*s

        @test NBodyMatrixElement(bcs, h, 1, 1) == ome(1, h, 1) + ome(2, h, 2) + ome(3, h, 3)
        @test NBodyMatrixElement(bcs, h, 1, 3) == -ome(3, h, 2)*s + ome(3, h, 5)
        @test NBodyMatrixElement(bcs, h, 1, 4) == ome(3, h, 6)*s

        @test NBodyMatrixElement(bcs, g, 1, 1) == (
            (ome([2,3], g, [2,3]) +
                ome([1,3], g, [1,3]) +
                ome([2,1], g, [2,1]))
            -
                (ome([2,3], g, [3,2]) +
                ome([1,3], g, [3,1]) +
                ome([2,1], g, [1,2])))

        @test NBodyMatrixElement(bcs, g, 1, 4) == (
            (ome([3,1], g, [6,1])*s +
                ome([2,3], g, [5,6]))
            -
                (ome([3,1], g, [1,6])*s +
                ome([2,3], g, [6,5])))

        @test NBodyMatrixElement(bcs, H, 1, 1) ==
            ome(1, h, 1) + ome(2, h, 2) + ome(3, h, 3) +
            ((ome([2,3], g, [2,3]) +
            ome([1,3], g, [1,3]) +
            ome([2,1], g, [2,1]))
             -
                 (ome([2,3], g, [3,2]) +
                 ome([1,3], g, [3,1]) +
                 ome([2,1], g, [1,2])))
    end

    @testset "Multi-configurational" begin
        s = ov(2,5)
        bcs = BitConfigurations([[1,2,3],[1,5,6]], [s])

        h = FieldFreeOneBodyHamiltonian()
        val = nothing
        err = @capture_err begin
            val = Matrix(bcs, h, verbosity=Inf)
        end
        @test val ==  [ome(1, h, 1) + ome(2, h, 2) + ome(3, h, 3) ome(3, h, 6)*s
                       ome(6, h, 3)*s' ome(1, h, 1) + ome(5, h, 5) + ome(6, h, 6)]
        @test occursin("[ Info: Generating energy expression\n", err)
    end
end
