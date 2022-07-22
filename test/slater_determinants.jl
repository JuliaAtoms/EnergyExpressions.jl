@testset "Slater determinants" begin
    c = sc"1s₀α 2s₀α 2p₋₁β"
    os = orbitals(c)
    s = SlaterDeterminant(c)
    @test length(s) == 3
    @test num_electrons(s) == 3
    @test orbitals(s) == os

    @test string(s) == "|1s₀α 2s₀α 2p₋₁β|"
    @test strdisplay(s) == "1s₀α(1)2s₀α(2)2p₋₁β(3) - 1s₀α(1)2s₀α(3)2p₋₁β(2) - 1s₀α(2)2s₀α(1)2p₋₁β(3) + 1s₀α(2)2s₀α(3)2p₋₁β(1) + 1s₀α(3)2s₀α(1)2p₋₁β(2) - 1s₀α(3)2s₀α(2)2p₋₁β(1)"

    @test strdisplay(SlaterDeterminant(collect(1:3))) ==
        "1(1)2(2)3(3) - 1(1)2(3)3(2) - 1(2)2(1)3(3) + 1(2)2(3)3(1) + 1(3)2(1)3(2) - 1(3)2(2)3(1)"
    @test strdisplay(SlaterDeterminant(collect(1:4))) == "|1 2 3 4|"

    as = s'
    @test as isa EnergyExpressions.AdjointSlaterDeterminant
    @test length(as) == 3
    @test string(as) == "[|1s₀α 2s₀α 2p₋₁β|]†"
end
