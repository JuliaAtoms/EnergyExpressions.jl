@testset "NBodyMatrixElements" begin
    @testset "NBodyTermFactor" begin
        @testset "OrbitalOverlap" begin
            s = ov(1,2)
            @test !iszero(s)
            @test s == ov(1,2)
            @test hash(s) == hash(ov(1,2))
            @test s' == ov(2,1)
            @test numbodies(s) == 0

            @test isdependent(s, 2)
            @test !isdependent(s, 1)

            @test !isdependent(s, Conjugate(2))
            @test isdependent(s, Conjugate(1))

            @test string(s) == "⟨1|2⟩"
        end

        @testset "OrbitalMatrixElement" begin
            Z = IdentityOperator{0}()
            A = IdentityOperator{1}()
            B = IdentityOperator{2}()

            @test_throws ArgumentError ome([1], A, [2, 3])
            @test_throws ArgumentError ome(1, B, 2)

            l = ome([], Z, [])
            m = ome(1, A, 2)
            n = ome([1,2], B, [2,3])

            @test !iszero(l)
            @test !iszero(m)
            @test !iszero(n)
            @test l == ome([], Z, [])
            @test m == ome(1, A, 2)
            @test n == ome([1,2], B, [2,3])
            @test hash(l) == hash(ome([], Z, []))
            @test hash(m) == hash(ome(1, A, 2))
            @test hash(n) == hash(ome([1,2], B, [2,3]))

            @test numbodies(l) == 0
            @test numbodies(m) == 1
            @test numbodies(n) == 2

            @test !isdependent(l, 1)

            @test isdependent(m, 2)
            @test !isdependent(m, 1)

            @test !isdependent(m, Conjugate(2))
            @test isdependent(m, Conjugate(1))

            @test string(l) == "⟨𝐈₀⟩"
            @test string(m) == "⟨1|𝐈₁|2⟩"
            @test string(n) == "⟨1 2|𝐈₂|2 3⟩"

            @test contract(n) == B
            @test contract(n, 1) == c(1, B, 2)
            @test contract(n, 2) == c(2, B, 3)
            @test contract(n, 1, 2) == c([1,2], B, [2,3])
        end
    end

    @testset "NBodyTerm" begin
        o = one(NBodyTerm)
        z = zero(NBodyTerm)

        @test o == one(o)
        @test isone(o)
        @test EnergyExpressions.isminusone(-o)

        @test z == zero(o)
        @test iszero(z)

        h = FieldFreeOneBodyHamiltonian()
        s = ov(1,2)

        m = ome(3,h,4)
        n = nbt(m, s)
        @test n == nbt(m, s)
        @test n == convert(NBodyTerm, m)*convert(NBodyTerm, s)
    end

    @testset "NBodyMatrixElement" begin
        o = one(NBodyMatrixElement)
        z = zero(NBodyMatrixElement)

        @test o == one(z)
        @test isone(o)

        @test z == zero(o)
        @test iszero(z)


        h = FieldFreeOneBodyHamiltonian()
        s = ov(1,2)
        m = ome(3,h,4)
        n = nbt(m, s)

        @test convert(NBodyMatrixElement, n) == nbme(n)
        @test convert(NBodyMatrixElement, 3) == nbme(3one(NBodyTerm))

        @test nbme(n) == n
        @test n == nbme(n)

        me = nbme(n)
        @test me ≈ (1+1e-9)*me
        @test me ≈ 1.001me atol=1e-3

        @test string(zero(NBodyMatrixElement)) == "0"
        @test string(0*(m+s)) == "0"
        @test string(2*(m+s)) == "2.0(3|4) + 2.0⟨1|2⟩"

        me = m + n + s
        me += me
        me += me
        me += me
        @test strdisplay(me, limit=true) == "(3|4) + (3|4)⟨1|2⟩ + ⟨1|2⟩ + (3|4) + (3|4)⟨1|2⟩ + ⟨1|2⟩ + … + (3|4) + (3|4)⟨1|2⟩ + ⟨1|2⟩ + (3|4) + (3|4)⟨1|2⟩ + ⟨1|2⟩"
    end

    @testset "Arithmetic" begin
        h = FieldFreeOneBodyHamiltonian()
        s = ov(1,2)

        @test 2s == nbt(s,f=2)
        @test s*2 == 2s
        @test -s == nbt(s,f=-1)
        @test -s == (-1)*s

        m = ome(3,h,4)
        n = nbt(m, s)

        @test 2n == nbt(m,s, f=2)
        @test n*2 == 2n

        @test n*(2n) == nbt(m,s,m,s, f=2) == nbt(s,m,s,m, f=2)
        @test -n == (-1)*n

        @test n*s == nbt(m,s,s, f=1)
        @test s*n == nbt(m,s,s, f=1)
        @test s*s == nbt(s,s, f=1)

        @test s + s == 2s
        @test s + m == nbme(s, m)
        @test 2s + m == nbme(2s, m)
        @test s + 2m == nbme(s, 2m)

        @test isdependent(n, 2)
        @test isdependent(n, 4)
        @test !isdependent(n, 3)
        @test !isdependent(n, 1)

        @test !isdependent(n, Conjugate(2))
        @test !isdependent(n, Conjugate(4))
        @test isdependent(n, Conjugate(3))
        @test isdependent(n, Conjugate(1))

        @test string(n) == "(3|4)⟨1|2⟩"
        @test string(-n) == "- (3|4)⟨1|2⟩"
        @test string(2n) == "2.0(3|4)⟨1|2⟩"
        @test string(-2n) == "- 2.0(3|4)⟨1|2⟩"
        @test strdisplay(n*n*n, limit=true) == "(3|4)⟨1|2⟩…(3|4)⟨1|2⟩"

        @test (m+s) + (m+s) == 2m + 2s
        @test (m+s) + s == m + 2s
        @test s + (m+s) == m + 2s
        @test (m+s) + 2s == m + 3s
        @test 2s + (m+s) == m + 3s

        @test m - s == nbme(m, -s)
        @test 2m - 2s == 2nbme(m, -s)
        @test 2*(m+s) == 2m + 2s
        @test (m+s)*2 == 2m + 2s
        @test_broken (m+s)*s == m*s + s*s
        @test (m+s)*s == s*m + s*s
        @test -(m+s) == -m - s
    end

    @testset "Matrix elements" begin
        # We test the implementation based on SlaterDeterminant:s
        # using the other implementation based on BitConfigurations.
        a = [1,2,6]
        b = [3,4,7]
        s = ov(6,7)

        bcs = BitConfigurations([a,b],[s])

        sa = SlaterDeterminant(a)
        sb = SlaterDeterminant(b)
        cfgs = [sa, sb]

        h = FieldFreeOneBodyHamiltonian()
        g = CoulombInteraction()

        H = h + g

        E1 = Matrix(bcs, H)
        E2 = nothing
        err = @capture_err begin
            E2 = Matrix(H, cfgs, [s], verbosity=Inf)
        end
        # The old implementation generates the opposite order of all
        # factors, and we do not yet know how to compare matrix
        # elements with flipped factor order, therefore we flip them
        # manually.
        for nv in E2.nzval
            for t in nv.terms
                reverse!(t.factors)
            end
        end

        @test E1 == E2
        @test occursin("[ Info: Generating energy expression\n", err)
    end
end
