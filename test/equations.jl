@testset "Equations" begin
    h = FieldFreeOneBodyHamiltonian()
    g = CoulombInteraction()

    z = zero(NBodyEquation)
    @test iszero(z)

    e = eq(:a, c(:a, g, :b))
    e2 = 2e
    @test !iszero(e)
    @test !iszero(e2)
    @test e2.factor.coeff == 2
    @test iszero(0e)

    @test -e == (-1)*e
    @test (-e).factor.coeff == -1

    @test string(e) == "+ 1[a|b]|a⟩"
    @test string(ceq(:a, c(:a, g, :b))) == "+ 1⟨a|[a|b]"

    ee = eq(:a, c(:a, g, :b)) + eq(:a, h) + eq(:b, c(:a, g, :a))

    @test -ee == -eq(:a, c(:a, g, :b)) - eq(:a, h) - eq(:b, c(:a, g, :a))

    @test iszero(zero(LinearCombinationEquation))
    @test iszero(zero(ee))
    @test iszero(ee + (-1)*ee)

    @test 0ee - ee == -ee
    @test ee - 0ee == ee

    @test ee - 2eq(:b, h) == eq(:a, c(:a, g, :b)) + eq(:a, h) + eq(:b, c(:a, g, :a)) - 2eq(:b, h)
    @test 2eq(:b, h) - ee == 2eq(:b, h) - eq(:a, c(:a, g, :b)) - eq(:a, h) - eq(:b, c(:a, g, :a))
    @test (eq(:a, c(:a, g, :b)) + eq(:a, h)) - eq(:a, c(:a, g, :b)) == eq(:a, h)
    @test (eq(:a, c(:a, g, :b)) + eq(:a, h)) - 2eq(:a, c(:a, g, :b)) == (-eq(:a, c(:a, g, :b)) + eq(:a, h))
    @test (eq(:a, c(:a, g, :b)) + eq(:a, h)) - eq(:b, h) == (eq(:a, c(:a, g, :b)) + eq(:a, h)) + (-1)*eq(:b, h)

    @test string(ee) == "+ 1[a|b]|a⟩ + 1ĥ₀|a⟩ + 1[a|a]|b⟩"

    @test strdisplay(eq(:a, c(:a, g, :b)) + eq(:a, h) + eq(:b, c(:a, g, :a)) +
        eq(:a, c(:b, g, :a)) + eq(:b, h) + eq(:b, c(:a, g, :b)),
                     limit=true) == "+ 1[a|b]|a⟩ + 1ĥ₀|a⟩ … + 1ĥ₀|b⟩ + 1[a|b]|b⟩"

    @test string(0ee) == "0"


    @testset "Comparison of equations" begin
        test_equations_equal(eq(:a, c(:a, g, :b)) + eq(:a, h) + eq(:b, c(:a, g, :a)),
                             eq(:a, h) + eq(:a, c(:a, g, :b)) + eq(:b, c(:a, g, :a)))

        test_equations_equal(eq(:a, h) + eq(:b, c(:a, g, :a)) + eq(:a, h),
                             2eq(:a, h) + eq(:b, c(:a, g, :a)))
    end

end
