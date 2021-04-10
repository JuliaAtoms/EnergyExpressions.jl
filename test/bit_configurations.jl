import EnergyExpressions: BitPattern, Orbitals, BitConfiguration,
    showbin, num_particles, orbital_numbers, get_orbitals,
    num_excitations, holes_mask, holes, particles_mask, particles,
    holes_particles_mask, holes_particles,
    relative_sign

using SparseArrays

@testset "Show binary values" begin
    @test string(BitPattern(1)) == "1"
    @test string(BitPattern(2)) == ".1"
    @test string(BitPattern(3)) == "11"
    @test string(BitPattern(5)) == "1.1"
    @test string(BitPattern(16)) == ".... 1"
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
            o = Orbitals(collect(1:10), [OrbitalOverlap(3,8)])
            @test o.overlaps isa SparseMatrixCSC
            @test size(o.overlaps) == (10,10)
            @test o.overlaps[3,8].factors[1] == OrbitalOverlap(3,8)
            @test o.overlaps[8,3].factors[1] == OrbitalOverlap(8,3)

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

        @test num_particles(a) == 2
        @test num_particles(b) == 3

        @test orbital_numbers(a) == [1,3]
        @test orbital_numbers(b) == [1,2,5]

        @test get_orbitals(a) == ["a", "c"]
        @test get_orbitals(b) == ["a", "b", "e"]

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

        @test relative_sign(b,d,holes_particles(b,d)...,verbosity=Inf) == 1
        @test relative_sign(b,e,holes_particles(b,e)...,verbosity=Inf) == -1
    end
end
