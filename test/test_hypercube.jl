using Test
using LinearAlgebra
using TightBindingLattice

@testset "hypercube" begin
    @testset "ExactLinearAlgebra" begin
        mat = Int[1 2 3;
                  4 5 2;
                  7 8 10]

        @testset "cofactor" begin
            results = Dict(
                (1,1) => [5 2; 8 10],
                (1,2) => [4 2; 7 10],
                (1,3) => [4 5; 7 8],
                (2,1) => [2 3; 8 10],
                (2,2) => [1 3; 7 10],
                (2,3) => [1 2; 7 8],
                (3,1) => [2 3; 5 2],
                (3,2) => [1 3; 4 2],
                (3,3) => [1 2; 4 5]
            )

            for i in 1:3, j in 1:3
                out = ones(Int, (2,2))
                ExactLinearAlgebra.get_cofactor_matrix_unsafe!(out, mat, i, j)
                @test out == results[i,j]
            end
        end

        @testset "determinant" begin

            @test_throws ArgumentError ExactLinearAlgebra.determinant([1 2 3; 4 5 6])
            @test_throws ArgumentError ExactLinearAlgebra.determinant(ones(Int, 0,0))
            @test ExactLinearAlgebra.determinant(7*ones(Int, 1, 1)) == 7

            @test abs(det(mat) - ExactLinearAlgebra.determinant(mat)) < sqrt(eps(Float64))
            for n in [5, 6, 7, 8]
                a = rand(-3:3, (n, n))
                det1 = det(a)
                det2 = ExactLinearAlgebra.determinant(a)
                @test abs(det1 - det2) < sqrt(eps(Float64))
            end
        end # testset determinant

        @testset "inverse" begin
            @test_throws ArgumentError ExactLinearAlgebra.inverse([1 2 3; 4 5 6])
            @test_throws ArgumentError ExactLinearAlgebra.inverse(ones(Int, 0, 0))
            @test ExactLinearAlgebra.inverse(7*ones(Int, 1, 1)) == ones(Int, 1, 1)//7

            @test maximum(abs.(inv(mat) - ExactLinearAlgebra.inverse(mat))) < sqrt(eps(Float64))
            for n in [5, 6, 7, 8]
                a = rand(-3:3, (n, n)) + 30*I
                inv1 = inv(a)
                inv2 = ExactLinearAlgebra.inverse(a)
                @test maximum(abs.(inv1 - inv2)) < sqrt(eps(Float64))
            end

        end
    end # testset ExactLinearAlgebra


    @testset "HypercubicLattice" begin
        @testset "failures" begin
            @test_throws DimensionMismatch HypercubicLattice([1 2 3; 4 5 6])
            @test_throws ArgumentError HypercubicLattice([1 1; 1 1])
            @test_throws ArgumentError HypercubicLattice([2 0; 0 2], [[0,0]])
            @test_throws ArgumentError HypercubicLattice([2 0; 0 2], [[0,0],[1,0],[0,1],[1,2]])
        end

        @testset "orthogonal" begin
            hypercube = HypercubicLattice([3 0; 0 3])
            hypercube2 = HypercubicLattice([3 0; 0 3],
                                                                          [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]])
            @test hypercube == hypercube2
            @test hypercube.inverse_scale_matrix == [1//3 0; 0 1//3]
            @test hypercube2.inverse_scale_matrix == [1//3 0; 0 1//3]

            for i in 0:6, j in 0:6
                @test hypercube.wrap([i, j]) == ([i÷3, j÷3], [i%3, j%3])
                @test hypercube2.wrap([i, j]) == ([i÷3, j÷3], [i%3, j%3])
            end
            let r = hcat([[i,j] for i in 0:6 for j in 0:6]...),
                    R = hcat([[i÷3,j÷3] for i in 0:6 for j in 0:6]...),
                    ρ = hcat([[i%3,j%3] for i in 0:6 for j in 0:6]...)
                    @test hypercube.wrap(r) == (R, ρ)
                    @test hypercube2.wrap(r) == (R, ρ)
            end

            @test dimension(hypercube) == 2
            @test isequiv(hypercube, HypercubicLattice([3 0; 3 3]))
        end

        @testset "nonorthogonal" begin
            #hypercube = HypercubicLattice([4 -2; 2 2])
            #hypercube = HypercubicLattice([3 0; 0 2])
            #@show hypercube.coordinates
            hypercube = HypercubicLattice([2 1; 1 2])
            @test hypercube.inverse_scale_matrix == [2//3 -1//3; -1//3 2//3]
            @test hypercube.coordinates == [[0,0], [2,2], [1,1]]
            @test dimension(hypercube) == 2
            @test hypercube.wrap([0,0]) == ([0,  0], [0, 0])
            @test hypercube.wrap([1,0]) == ([0, -1], [2, 2])
            @test hypercube.wrap([2,0]) == ([1, -1], [1, 1])
            @test hypercube.wrap([3,0]) == ([2, -1], [0, 0])
        end

        # @testset "nonorthogonal-gen" begin
        #     hypercube = HypercubicLattice([3 0; 0 2])
        #     @test hypercube.inverse_scale_matrix == [1//3 0; 0 1//2]
        #     hypercube2 = orthogonalize(hypercube)
        #     @test hypercube.coordinates == [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1]]
        #     @test hypercube2.coordinates == [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1]]   # generator: [1,1]
        #     for i in -3:3, j in -3:3
        #         @test hypercube.wrap([i,j]) == hypercube2.wrap([i,j])
        #     end
        # end
    end
end # testset enlargement
