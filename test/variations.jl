@testset "Calculus-of-variations" begin
    h = FieldFreeOneBodyHamiltonian()
    g = CoulombInteraction()
    I₁ = IdentityOperator{1}()
    I₃ = IdentityOperator{3}()
    I₄ = IdentityOperator{4}()

    @testset "Comparison of equations" begin
        test_equations_equal(eq(:a, c(:a, g, :b)) + eq(:a, h) + eq(:b, c(:a, g, :a)),
                             eq(:a, h) + eq(:a, c(:a, g, :b)) + eq(:b, c(:a, g, :a)))

        test_equations_equal(eq(:a, h) + eq(:b, c(:a, g, :a)) + eq(:a, h),
                             2eq(:a, h) + eq(:b, c(:a, g, :a)))
    end

    @testset "Various variations" begin
        extra = EnergyExpressions.OrbitalMatrixElement([:a,:a], g, [:a,:b])
        ex_extra = eq(:a, c(:a, g, :b)) + eq(:b, c(:a, g, :a))
        cex_extra = ceq(:a, c(:a, g, :b))

        for (a,op,b,ovs,ex,cex) in [
            ([:a], h, [:a], (), eq(:a, h), ceq(:a, h)),
            ([:a], h, [:b], (), eq(:b, h), 0),
            ([:b], h, [:a], (), 0, ceq(:b, h)),
            ([:a], h, [:a], ((:a,:a),),
             eq(:a, h, nbt(ov(:a, :a))) + eq(:a, I₁, nbt(ome(:a, h, :a))),
             ceq(:a, h, nbt(ov(:a, :a))) + ceq(:a, I₁, nbt(ome(:a, h, :a)))),
            ([:a], h, [:b], ((:a,:a),),
             eq(:b, h, nbt(ov(:a, :a))) + eq(:a, I₁, nbt(ome(:a, h, :b))),
             ceq(:a, I₁, nbt(ome(:a, h, :b)))),
            ([:b], h, [:a], ((:a,:a),),
             eq(:a, I₁, nbt(ome(:b, h, :a))),
             ceq(:b, h, nbt(ov(:a, :a))) + ceq(:a, I₁, nbt(ome(:b, h, :a)))),
            ([:a], h, [:a], ((:a,:b),),
             eq(:a, h, nbt(ov(:a, :b))) + eq(:b, I₁, nbt(ome(:a, h, :a))),
             ceq(:a, h, nbt(ov(:a, :b)))),
            ([:a], h, [:b], ((:a,:b),),
             eq(:b, h, nbt(ov(:a, :b))) + eq(:b, I₁, nbt(ome(:a, h, :b))),
             0),
            ([:b], h, [:a], ((:a,:b),),
             eq(:b, I₁, nbt(ome(:b, h, :a))),
             ceq(:b, h, nbt(ov(:a, :b)))),
            ([:a], h, [:a], ((:c,:d),),
             eq(:a, h, nbt(ov(:c, :d))),
             ceq(:a, h, nbt(ov(:c, :d)))),
            ([:a], h, [:b], ((:c,:d),),
             eq(:b, h, nbt(ov(:c, :d))),
             0),
            ([:b], h, [:a], ((:c,:d),),
             0,
             ceq(:b, h, nbt(ov(:c, :d)))),
            ([:a,:a], g, [:a,:a], (),
             2eq(:a, c(:a, g, :a)),
             2ceq(:a, c(:a, g, :a))),
            ([:a,:b], g, [:a,:b], (),
             eq(:a, c(:b, g, :b)),
             ceq(:a, c(:b, g, :b))),
            ([:a,:b], g, [:b,:a], (),
             eq(:b, c(:b, g, :a)),
             ceq(:b, c(:a, g, :b))),
            ([:a,:a], g, [:a,:b], (),
             eq(:a, c(:a, g, :b)) + eq(:b, c(:a, g, :a)),
             ceq(:a, c(:a, g, :b))),
            ([:a,:a], g, [:a,:b], ((:a,:a),),
             (eq(:a, c(:a, g, :b), nbt(ov(:a, :a))) + eq(:b, c(:a, g, :a), nbt(ov(:a, :a))) +
                 eq(:a, I₁, nbt(ome([:a,:a], g, [:a,:b])))),
             (ceq(:a, c(:a, g, :b), nbt(ov(:a, :a))) +
                 ceq(:a, I₁, nbt(ome([:a,:a], g, [:a,:b]))))),
            ([:a,:a], g, [:a,:b], ((:c,:d),),
             (eq(:a, c(:a, g, :b), nbt(ov(:c, :d))) + eq(:b, c(:a, g, :a), nbt(ov(:c, :d)))),
             (ceq(:a, c(:a, g, :b), nbt(ov(:c, :d))))),
            ([:a,:a,:a], I₃, [:a,:a,:a], (),
             3eq(:a, c([:a,:a], I₃, [:a,:a])),
             3ceq(:a, c([:a,:a], I₃, [:a,:a]))),
            ([:a,:a,:a], I₃, [:a,:b,:a], (),
             eq(:a, c([:a,:a], I₃, [:b,:a])) + eq(:b, c([:a,:a], I₃, [:a,:a])) + eq(:a, c([:a,:a], I₃, [:a,:b])),
             ceq(:a, c([:a,:a], I₃, [:b,:a])) + ceq(:a, c([:a,:a], I₃, [:a,:b]))),
            ([:a,:a,:a], I₃, [:a,:b,:c], (),
             eq(:a, c([:a,:a], I₃, [:b,:c])) + eq(:b, c([:a,:a], I₃, [:a,:c])) + eq(:c, c([:a,:a], I₃, [:a,:b])),
             ceq(:a, c([:a,:a], I₃, [:b,:c]))),
            ([:a,:a,:a], I₃, [:a,:a,:a], ((:a,:a),),
             3eq(:a, c([:a,:a], I₃, [:a,:a]), nbt(ov(:a, :a))) + eq(:a, I₁, nbt(ome([:a, :a, :a], I₃, [:a, :a, :a]))),
             3ceq(:a, c([:a,:a], I₃, [:a,:a]), nbt(ov(:a, :a))) + ceq(:a, I₁, nbt(ome([:a, :a, :a], I₃, [:a, :a, :a])))),
            ([:a, :a, :a, :a], I₄, [:a, :a, :b, :a], (),
             2eq(:a, c([:a, :a, :a], I₄, [:a, :b, :a])) + eq(:b, c([:a, :a, :a], I₄, [:a, :a, :a])) + eq(:a, c([:a, :a, :a], I₄, [:a, :a, :b])),
             2ceq(:a, c([:a, :a, :a], I₄, [:a, :b, :a])) + ceq(:a, c([:a, :a, :a], I₄, [:a, :a, :b]))),
            ([:a, :a, :a, :a], I₄, [:a, :b, :a, :a], (),
             eq(:a, c([:a, :a, :a], I₄, [:b, :a, :a])) + eq(:b, c([:a, :a, :a], I₄, [:a, :a, :a])) + 2eq(:a, c([:a, :a, :a], I₄, [:a, :b, :a])),
             ceq(:a, c([:a, :a, :a], I₄, [:b, :a, :a])) + 2ceq(:a, c([:a, :a, :a], I₄, [:a, :b, :a])))
            ]
            local ome = EnergyExpressions.OrbitalMatrixElement(a, op, b)
            local nbt = isempty(ovs) ? ome : EnergyExpressions.NBodyTerm(vcat(ome, [OrbitalOverlap(ov...) for ov in ovs]), 1)
            nbme1 = EnergyExpressions.NBodyMatrixElement([nbt])
            nbme2 = EnergyExpressions.NBodyMatrixElement([nbt, extra])

            ex = iszero(ex) ? zero(NBodyEquation) : ex
            cex = iszero(cex) ? zero(NBodyEquation) : cex

            test_variations_equal(nbme1, Conjugate(:a), ex)
            test_variations_equal(nbme1, :a, cex)

            test_variations_equal(nbme2, Conjugate(:a), ex + ex_extra)
            test_variations_equal(nbme2, :a, cex + cex_extra)
        end
    end
end
