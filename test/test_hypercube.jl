using Test
using LinearAlgebra
using TightBindingLattice

@testset "hypercube" begin


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
