@testset "Loop macros" begin
    @testset "@above_diagonal_loop" begin
        v = Vector{NTuple{3,Int}}()
        EnergyExpressions.@above_diagonal_loop 3 i 5 begin
            @test i_1 > i_2 > i_3
            push!(v, (i_1, i_2, i_3))
        end
        @test length(v) == 10
        @test v[1] == (3,2,1)
        @test v[end] == (5,4,3)
    end

    @testset "@anti_diagonal_loop" begin
        v = Vector{NTuple{4,Int}}()
        EnergyExpressions.@anti_diagonal_loop 4 i 10 begin
            push!(v, Base.Cartesian.@ntuple(4, i))
        end
        @test length(v) == 5040
        @test allunique(v)
        @test all(e -> length(e) == 4, v)
        @test all(allunique, v)
        @test all(e -> all(in(1:10), extrema(e)), v)
    end
end
