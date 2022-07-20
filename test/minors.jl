# Used to generate reference data to test the linear complexity
# algorithm against; the naïve algorithm is correct, but is of
# factorial complexity and cannot be used for more than ~20 orbitals.
@generated function naive_nonzero_minors(op::NBodyOperator{N}, overlap, numorbitals) where N
    quote
        ks = Vector{Vector{Int}}()
        ls = Vector{Vector{Int}}()

        @above_diagonal_loop $N k numorbitals begin
            kN = zeros(Int, $N)
            Base.Cartesian.@nexprs $N d -> kN[d] = k_d

            @above_diagonal_loop $N l numorbitals begin
                lN = zeros(Int, $N)
                Base.Cartesian.@nexprs $N d -> lN[d] = l_d

                Dkl = cofactor(kN, lN, overlap) # Determinantal overlap of all other orbitals
                iszero(Dkl) && continue
                push!(ks, kN)
                push!(ls, lN)
            end
        end
        ks, ls
    end
end

struct ConcreteNBodyOperator{N} <: NBodyOperator{N} end

function test_find_minors(cfga, cfgb, overlaps=OrbitalOverlap[])
    a = SlaterDeterminant(cfga)
    b = SlaterDeterminant(cfgb)

    S = overlap_matrix(a, b, overlaps)

    h = ConcreteNBodyOperator{1}()
    g = ConcreteNBodyOperator{2}()

    errors = 0
    for (op,label) in [(h,"γ"), (g,"Γ")]
        ks,ls = naive_nonzero_minors(op, S, length(a))
        N = numbodies(op)
        ks′,ls′ = nonzero_minors(N, S)
        kls′ = collect(zip(sort.(ks′),sort.(ls′)))
        for (k,l) in zip(sort.(ks),sort.(ls))
            if (k,l) ∉ kls′
                errors += 1
                kl = nice_string(k)
                ll = nice_string(l)
                println("$label($kl|$ll) missing from new algorithm")
            end
        end
    end
    @test errors == 0
end

function find_diagonal_blocks_test_case(A, blocks)
    @test find_diagonal_blocks(A) == blocks
    @test find_diagonal_blocks(SparseMatrixCSC(A')) == blocks
end

function test_cofactors(N, ii, jj)
    orbitals = 1:N
    excite = i -> i+100
    virtuals = map(excite, orbitals)

    configurations = [SlaterDeterminant(vcat(1:i-1,excite(i),i+1:N))
                      for i=1:N]
    overlaps = map(((a,b),) -> OrbitalOverlap(a,b),
                   reduce(vcat, [collect(combinations(virtuals,2)),
                                 [[a,a] for a in virtuals]]))

    results = Vector{EnergyExpressions.NBodyMatrixElement}()

    for (i,j) = [(ii,jj),(jj,ii)]
        a = configurations[i]
        b = configurations[j]
        S = overlap_matrix(a, b, overlaps)
        ks,ls = nonzero_minors(1, S)
        for (k,l) in zip(ks,ls)
            push!(results, cofactor(k, l, S))
        end
    end
    results
end

@testset "Minor determinants" begin
    @testset "Find diagonal blocks" begin
        find_diagonal_blocks_test_case((sparse([3,1,2], 1:3, ones(Int, 3))), [1:3])
        find_diagonal_blocks_test_case((sparse([1,4,2,3], 1:4, ones(Int, 4))), [1:1, 2:4])
        find_diagonal_blocks_test_case((sparse([1,2,5,3,4], 1:5, ones(Int, 5))), [1:1, 2:2, 3:5])
        find_diagonal_blocks_test_case((sparse([1,2,3,6,4,5], 1:6, ones(Int, 6))), [1:1, 2:2, 3:3, 4:6])
    end

    @testset "Cofactors" begin
        @test test_cofactors(7,4,7) == [NBodyMatrixElement([-OrbitalOverlap(104,107)]),
                                          NBodyMatrixElement([-OrbitalOverlap(107,104)])]
        @test test_cofactors(8,5,8) == [NBodyMatrixElement([-OrbitalOverlap(105,108)]),
                                          NBodyMatrixElement([-OrbitalOverlap(108,105)])]
        @test test_cofactors(9,6,9) == [NBodyMatrixElement([-OrbitalOverlap(106,109)]),
                                          NBodyMatrixElement([-OrbitalOverlap(109,106)])]
        @test test_cofactors(9,7,9) == [NBodyMatrixElement([-OrbitalOverlap(107,109)]),
                                          NBodyMatrixElement([-OrbitalOverlap(109,107)])]
        @test test_cofactors(54,45,54) == [NBodyMatrixElement([-OrbitalOverlap(145,154)]),
                                             NBodyMatrixElement([-OrbitalOverlap(154,145)])]
    end

    @testset "Base cases" begin
        test_find_minors([:a], [:a])
        test_find_minors([:a,:b], [:a,:b])
        test_find_minors([:a,:b], [:b,:a])
        test_find_minors([:a,:b], [:a,:c])
        test_find_minors([:a,:b], [:c,:a])
        test_find_minors([:a,:b,:c], [:a,:b,:c])
        test_find_minors([:a,:b,:c], [:a,:b,:d])
        test_find_minors([:a,:b,:c], [:a,:b,:d], [OrbitalOverlap(:a,:b)])
        test_find_minors([:a,:b,:c], [:a,:b,:d], [OrbitalOverlap(:a,:c)])
        test_find_minors([:a,:b,:c], [:a,:b,:d], [OrbitalOverlap(:b,:d),OrbitalOverlap(:c,:d)])
        test_find_minors([:a,:b,:c], [:a,:b,:d], [OrbitalOverlap(:a,:c),OrbitalOverlap(:c,:d)])
        test_find_minors([:a,:b,:c], [:a,:b,:d], [OrbitalOverlap(:a,:d),OrbitalOverlap(:c,:d)])
        test_find_minors([:a,:b,:c], [:b,:c,:a])
        test_find_minors([:a,:b,:c], [:c,:b,:a])
        test_find_minors([:a,:b,:c], [:b,:c,:a], [OrbitalOverlap(:a,:c)])
        test_find_minors([:a,:b,:c], [:a,:d,:e])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:c,:d])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:c,:d], [OrbitalOverlap(:c,:d)])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:c,:e], [OrbitalOverlap(:c,:d)])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:c,:e], [OrbitalOverlap(:c,:e)])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:e,:f], [OrbitalOverlap(:c,:b),OrbitalOverlap(:d,:b)])
        test_find_minors([:a,:b,:c,:d], [:a,:b,:e,:f], [OrbitalOverlap(:a,:e),OrbitalOverlap(:c,:f),OrbitalOverlap(:d,:b)])
        test_find_minors([:b,:c,:d,:k_a], [:a,:c,:d,:k_b], [OrbitalOverlap(:k_a,:k_b)])
        test_find_minors([:k_a,:b,:c,:d], [:a,:k_b,:c,:d], [OrbitalOverlap(:k_a,:k_b)])
        test_find_minors([:a, :b, :c, :d], [:a, :e, :f, :g])
        test_find_minors([:a,:b,:c,:d], [:d,:c,:b,:a])
    end

    @testset "Atomic cases" begin
        h = ConcreteNBodyOperator{1}()
        g = ConcreteNBodyOperator{2}()

        for (cfga,cfgb,exp1,exp2) in [
            (c"1s2",c"1s2",2,1),
            (c"1s2",c"1s 2p",1,1),
            (c"1s2 2s2",c"1s2 2s2",4,binomial(4,2)),
            (c"1s2 2s2 2p6",c"1s2 2s2 2p6",10,binomial(10,2)),
            (c"1s2 2s2 2p6 3s2",c"1s2 2s2 2p6 3s2",12,binomial(12,2)),
            (c"1s2 2s2 2p6 3s2 3p6",c"1s2 2s2 2p6 3s2 3p6",18,binomial(18,2)),
            (c"[Ar]",c"[Ar]",18,binomial(18,2)),
            (c"[Kr]",c"[Kr]",36,binomial(36,2)),
            (c"[Xe]",c"[Xe]",54,binomial(54,2))
        ]
            a = SlaterDeterminant(spin_configurations(cfga)[1])
            b = SlaterDeterminant(spin_configurations(cfgb)[1])

            S = overlap_matrix(a, b)

            for (op,ex) in [(h,exp1),(g,exp2)]
                ks,ls = nonzero_minors(numbodies(op), S)
                @test length(ks) == length(ls) == ex
            end
        end
    end
end
