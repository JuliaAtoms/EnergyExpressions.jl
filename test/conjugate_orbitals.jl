@testset "Conjugated orbitals" begin
    for (o,label) in ((o"1s", "1s†"), (so"1s(0,α)", "1s₀α†"), (rso"3d-(-1/2)", "3d-(-1/2)†"))
        co = conj(o)
        @test co == Conjugate(o)
        @test conj(co) == o
        @test symmetry(co) == symmetry(o)
        @test string(co) == label
    end
end
