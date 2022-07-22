@testset "KeyTracker" begin
    kt = EnergyExpressions.KeyTracker{Int}()

    for (i,j) in enumerate(10:-1:1)
        @test get!(kt, j) == i
    end

    @test keys(kt) == 10:-1:1
end
