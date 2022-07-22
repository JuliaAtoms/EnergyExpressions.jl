@testset "Invariant sets" begin
    @test invariant_sets(sparse(Matrix(I, 3,3))) == [[1], [2], [3]]
    @test invariant_sets(sparse([1 1 0; 1 1 0; 0 0 1])) == [[1, 2], [3]]
    # Unsure if we do not want state 3 to be its own set
    @test invariant_sets(sparse([0 1 0; 1 0 0; 0 0 0])) == [[1, 2]]
end
