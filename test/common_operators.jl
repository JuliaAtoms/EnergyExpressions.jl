import EnergyExpressions: OrbitalMatrixElement

@testset "Common operators" begin
    @testset "OneBodyHamiltonian" begin
        h = OneBodyHamiltonian()
        @test string(h) == "ĥ"
        ahb = OrbitalMatrixElement(["a"],h,["b"])
        @test string(ahb) == "(a|b)"
        @test !iszero(ahb)

        @test !iszero(OrbitalMatrixElement([so"1s₀α"], h, [so"2s₀α"]))
        @test !iszero(OrbitalMatrixElement([so"1s₀α"], h, [so"2p₀α"]))
        @test iszero(OrbitalMatrixElement([so"1s₀α"], h, [so"2p₀β"]))
    end

    @testset "FieldFreeOneBodyHamiltonian" begin
        h₀ = FieldFreeOneBodyHamiltonian()
        @test string(h₀) == "ĥ₀"
        ah₀b = OrbitalMatrixElement(["a"],h₀,["b"])
        @test string(ah₀b) == "(a|b)"
        @test !iszero(ah₀b)

        @test !iszero(OrbitalMatrixElement([so"1s₀α"], h₀, [so"2s₀α"]))
        @test iszero(OrbitalMatrixElement([so"1s₀α"], h₀, [so"2p₀α"]))
        @test iszero(OrbitalMatrixElement([so"1s₀α"], h₀, [so"2p₀β"]))
    end

    @testset "CoulombInteraction" begin
        g = CoulombInteraction()
        @test hash(g) == hash(nothing)
        @test string(g) == "ĝ"

        abcd = OrbitalMatrixElement(["a","b"],g,["c","d"])
        @test string(abcd) == "[a b|c d]"

        abab = OrbitalMatrixElement(["a","b"],g,["a","b"])
        @test string(abab) == "F(a,b)"

        abba = OrbitalMatrixElement(["a","b"],g,["b","a"])
        @test string(abba) == "G(a,b)"

        @test !iszero(OrbitalMatrixElement([so"1s₀α",so"2p₁β"],g,[so"1s₀α",so"2p₁β"]))
        @test iszero(OrbitalMatrixElement([so"1s₀α",so"2p₁β"],g,[so"2p₁β",so"1s₀α"]))

        bd = ContractedOperator(NTuple{1,String}(("b",)),g,NTuple{1,String}(("d",)))
        @test bd isa CoulombPotential
        @test string(bd) == "[b|d]"
    end
end
