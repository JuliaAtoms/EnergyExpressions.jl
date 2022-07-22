@testset "Multi-configurational equations" begin
    @testset "pushifmissing!" begin
        v = [1,2,3,4]
        @test pushifmissing!(v, 3) == 3
        @test length(v) == 4

        @test pushifmissing!(v, 5) == 5
        @test length(v) == 5
    end

    @testset "One-body operator" begin
        gst = collect(1:3)
        cfgs = [vcat(gst[1:i-1],100+i,gst[i+1:end])
                for i in eachindex(gst)]
        bcs = BitConfigurations(vcat([gst], cfgs))

        h = FieldFreeOneBodyHamiltonian()
        E = Matrix(bcs, h)

        os = orbitals(bcs)
        norb = length(os)

        eqs = nothing
        err = @capture_err begin
            eqs = diff(E, Conjugate.(os), verbosity=Inf)
        end

        @test length(eqs.equations) == norb
        @test isempty(eqs.integrals)

        @test occursin("OrbitalEquation(1): ", string(eqs.equations[1]))

        term_refs = [[(1,1,1,1), (1,2,1,101), (3,3,1,1), (4,4,1,1)],
                     [(1,1,1,2), (2,2,1,2), (1,3,-1,102), (4,4,1,2)],
                     [(1,1,1,3), (2,2,1,3), (3,3,1,3), (1,4,1,103)],
                     [(2,1,1,1), (2,2,1,101)],
                     [(3,1,-1,2), (3,3,1,102)],
                     [(4,1,1,3), (4,4,1,103)]]

        for (i,e) in enumerate(eqs.equations)
            @test e.orbital == os[i]
            @test length(e.terms) == 1
            integral,ts = first(e.terms)
            @test integral == 0
            term_ref = [EnergyExpressions.MCTerm(i,j,coeff,h,so,Int[])
                        for (i,j,coeff,so) in term_refs[i]]
            @test ts == term_ref
        end

        @test occursin("[ Info: Deriving equations for $(norb) orbitals\n", err)
    end

    @testset "Two-body operator" begin
        cfgs = [[1,2],[1,3]]
        bcs = BitConfigurations(cfgs)

        g = CoulombInteraction()
        E = Matrix(bcs, g)

        os = [2,3]
        norb = length(os)

        eqs = diff(E, Conjugate.(os))

        @test length(eqs.equations) == norb
        @test length(eqs.integrals) == 3
        @test c(1, g, 1) ∈ eqs.integrals
        i11 = findfirst(isequal(c(1, g, 1)), eqs.integrals)
        @test c(1, g, 2) ∈ eqs.integrals
        i12 = findfirst(isequal(c(1, g, 2)), eqs.integrals)
        @test c(1, g, 3) ∈ eqs.integrals
        i13 = findfirst(isequal(c(1, g, 3)), eqs.integrals)

        e1 = eqs.equations[1]
        e2 = eqs.equations[2]

        @test e1.orbital == 2
        @test e2.orbital == 3

        @test occursin("OrbitalEquation(2): ", string(e1))
        @test occursin("OrbitalEquation(3): ", string(e2))

        z = zero(LinearCombinationEquation)
        @test e1.equation == [-eq(1, c(1,g,2)) + eq(2, c(1,g,1)) -eq(1, c(1, g, 3)) + eq(3, c(1, g, 1)); z z]
        @test e2.equation == [z z; -eq(1, c(1,g,2)) + eq(2, c(1,g,1)) -eq(1, c(1, g, 3)) + eq(3, c(1, g, 1))]

        @test (0 => [EnergyExpressions.MCTerm(1,1,-1,c(1,g,2),1,Int[])]) ∈ e1.terms
        @test (i11 => [EnergyExpressions.MCTerm(1,1,1,c(1,g,1),2,Int[]),
                       EnergyExpressions.MCTerm(1,2,1,c(1,g,1),3,Int[])]) ∈ e1.terms
        @test (i13 => [EnergyExpressions.MCTerm(1,2,-1,c(1,g,3),1,Int[])]) ∈ e1.terms

        @test (0 => [EnergyExpressions.MCTerm(2,2,-1,c(1,g,3),1,Int[])]) ∈ e2.terms
        @test (i11 => [EnergyExpressions.MCTerm(2,1,1,c(1,g,1),2,Int[]),
                       EnergyExpressions.MCTerm(2,2,1,c(1,g,1),3,Int[])]) ∈ e2.terms
        @test (i12 => [EnergyExpressions.MCTerm(2,1,-1,c(1,g,2),1,Int[])]) ∈ e2.terms
    end

    @testset "Modify energy expression as-you-go" begin
        g = CoulombInteraction()
        E = sparse(reshape(NBodyMatrixElement[1ome([1,2],g,[3,4])],1,1))

        ref1 = reshape(LinearCombinationEquation[eq(3, c(2, g, 4))],1,1)
        ref2 = reshape(LinearCombinationEquation[eq(4, c(1, g, 3))],1,1)

        os = 1:2
        norb = length(os)

        eqs = diff(E, Conjugate.(os))
        @test eqs.integrals == [c(2, g, 4), c(1, g, 3)]
        @test length(eqs.equations) == 2
        @test eqs.equations[1].orbital == 1
        @test eqs.equations[2].orbital == 2
        @test eqs.equations[1].equation == ref1
        @test eqs.equations[2].equation == ref2
        @test eqs.equations[1].terms == [1 => [EnergyExpressions.MCTerm(1, 1, 1, c(2, g, 4), 3, Int[])]]
        @test eqs.equations[2].terms == [2 => [EnergyExpressions.MCTerm(1, 1, 1, c(1, g, 3), 4, Int[])]]

        # Now we instead modify the energy expression during the
        # calculus of variations, to avoid double-counting the
        # contribution of the two-electron integral, if we use the
        # orbital equations to compute the total energy.
        eqs2 = nothing
        err = @capture_err begin
            eqs2 = diff(E, Conjugate.(os), verbosity=Inf) do E,corb
                rows = rowvals(E)
                vals = nonzeros(E)
                m,n = size(E)
                for j = 1:n
                    for k in nzrange(E, j)
                        i = rows[k]
                        v = vals[k]
                        vals[k] = NBodyMatrixElement(filter(t -> !isdependent(t, corb), v.terms))
                    end
                end
            end
        end

        @test eqs2.integrals == [c(2, g, 4)]
        @test length(eqs2.equations) == 2
        @test eqs2.equations[1] == eqs.equations[1]
        @test iszero(eqs2.equations[2].equation)
        @test isempty(eqs2.equations[2].terms)

        @test occursin("[ Info: Deriving equations for $(norb) orbitals\n", err)
    end
end
