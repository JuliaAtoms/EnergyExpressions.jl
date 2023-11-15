@testset "LockedDict" begin
    ld = EnergyExpressions.LockedDict{Int,Int}()
    @test isempty(ld)

    Threads.@threads for i = 1:1_000
        n = mod(i, 10)
        @test get!(ld, n, () -> n) == n
    end

    for (k,v) in ld
        @test k == v
    end

    @test length(ld) == 10
    @test sort(collect(pairs(ld))) == [i => i for i = 0:9]
end
