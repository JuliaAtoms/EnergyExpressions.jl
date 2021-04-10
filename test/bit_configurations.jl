import EnergyExpressions: BitPattern, Orbitals
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
end
