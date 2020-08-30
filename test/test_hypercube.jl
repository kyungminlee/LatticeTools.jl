using Test
using LatticeTools


@testset "Hypercube" begin

    @testset "failures" begin
        @test_throws DimensionMismatch Hypercube([1 2 3; 4 5 6])
        @test_throws ArgumentError Hypercube([1 1; 1 1])
    end

    @testset "orthogonal" begin
        hypercube = Hypercube([3 0; 0 3])
        @test dimension(hypercube) == 2
        @test isequiv(hypercube, Hypercube([3 0; 3 3]))
        @test !isequiv(hypercube, Hypercube([3 0; 1 3]))
        @test !isequiv(hypercube, Hypercube([6 0; 0 6]))
        @test !isequiv(Hypercube([6 0; 0 6]), hypercube)

        for i in 0:6, j in 0:6
            @test hypercube.wrap([i, j]) == ([i÷3, j÷3], [i%3, j%3])
        end
        let r = hcat([[i,j] for i in 0:6 for j in 0:6]...),
            R = hcat([[i÷3,j÷3] for i in 0:6 for j in 0:6]...),
            ρ = hcat([[i%3,j%3] for i in 0:6 for j in 0:6]...)
            @test hypercube.wrap(r) == (R, ρ)
        end
    end

    @testset "nonorthogonal" begin
        hypercube = Hypercube([2 1; 1 2])
        @test hypercube.inverse_shape_matrix == [2//3 -1//3; -1//3 2//3]
        @test dimension(hypercube) == 2
        @test hypercube.wrap([0,0]) == ([0,  0], [0, 0])
        @test hypercube.wrap([1,0]) == ([0, -1], [2, 2])
        @test hypercube.wrap([2,0]) == ([1, -1], [1, 1])
        @test hypercube.wrap([3,0]) == ([2, -1], [0, 0])
    end

    @testset "generator" begin
        cube = Hypercube([4 0; 0 3])
        generate_coordinates(cube, [ 1 4; -1 -3])
        generate_coordinates(cube, [ 1 0; 0 1])
        @test_throws DimensionMismatch generate_coordinates(cube, [ 1 0 0; 0 1 0]) # not square
        @test_throws ArgumentError generate_coordinates(cube, [ 2 0; 0 1]) # unimodular
        @test_throws ArgumentError generate_coordinates(cube, [ 1 0; 1 1])
        # c2 = generate_coordinates(cube, [1 0; 0 1])
    end

end
