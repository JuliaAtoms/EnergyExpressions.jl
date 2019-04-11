import EnergyExpressions: @above_diagonal_loop, cofactor, nonzero_minors

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

@testset "Minor determinants" begin
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
