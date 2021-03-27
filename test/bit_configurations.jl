import EnergyExpressions: BitPattern

@testset "Show binary values" begin
    @test string(BitPattern(1)) == "1"
    @test string(BitPattern(2)) == ".1"
    @test string(BitPattern(3)) == "11"
    @test string(BitPattern(5)) == "1.1"
    @test string(BitPattern(16)) == ".... 1"
end
