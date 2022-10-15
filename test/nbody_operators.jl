@testset "N-body operators" begin
    I₁ = IdentityOperator{1}()
    @test ishermitian(I₁)
    @test I₁ isa OneBodyOperator
    @test string(I₁) == "𝐈₁"

    I₂ = IdentityOperator{2}()
    @test I₂ isa TwoBodyOperator

    @test adjoint(I₁) == I₁

    h = FieldFreeOneBodyHamiltonian()
    @test h' == AdjointOperator(h)
    @test (h')' == h'' == h
    @test string(h') == "ĥ₀†"

    @test hash(h') == hash(AdjointOperator(h))
    @test numbodies(h') == 1

    @testset "LinearCombinationOperator" begin
        L = I₁ + I₂
        @test !iszero(L)
        @test iszero(zero(L))
        @test iszero(zero(typeof(L)))

        @test L.operators == [I₁ => 1, I₂ => 1]
        @test (I₁ - I₂).operators == [I₁ => 1, I₂ => -1]
        @test (2I₁ + I₂).operators == [I₁ => 2, I₂ => 1]
        @test (2I₁ - I₂).operators == [I₁ => 2, I₂ => -1]
        @test (I₁ + 2I₂).operators == [I₁ => 1, I₂ => 2]
        @test (I₁ - im*I₂).operators == [I₁ => 1, I₂ => -im]
        @test (2I₁ + 2I₂).operators == [I₁ => 2, I₂ => 2]
        @test (2I₁ - 2I₂).operators == [I₁ => 2, I₂ => -2]
        @test (I₁*2).operators == [I₁ => 2]
        @test (2*(2I₁)).operators == [I₁ => 4]
        @test ((2I₁)*2).operators == [I₁ => 4]

        @test L == I₁ + I₂
        @test hash(L) == hash(I₁ + I₂)
        @test numbodies(L) == 2

        @test string(L) == "𝐈₁ + 𝐈₂"
        @test string(I₁ + 2I₂) == "𝐈₁ + 2.0𝐈₂"
        @test string(I₁ - im*I₂) == "𝐈₁ - im𝐈₂"
        @test string(zero(L)) == "0"

        @test filter(Base.Fix2(isa, OneBodyOperator), L).operators == [I₁ => 1]

        @test -I₁ == (-1)*I₁
        @test -I₂ == (-1)*I₂
        @test -L == -I₁ - I₂

        @test adjoint(L) == L
    end

    @testset "ContractedOperator" begin
        co = ContractedOperator((:a,), I₁, (:b,))

        @test co == ContractedOperator((:a,), I₁, (:b,))

        @test :b ∈ co
        @test :a ∉ co

        @test Conjugate(:b) ∉ co
        @test Conjugate(:a) ∈ co

        @test_throws ArgumentError ContractedOperator((:a,:c), I₁, (:b,:d))

        @test string(co) == "⟨a|𝐈₁|b⟩"
        @test string(ContractedOperator((:a,:c), I₂, (:b,:d))) == "⟨a c|𝐈₂|b d⟩"
    end
end
