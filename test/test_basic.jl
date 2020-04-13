using Test
using LinearAlgebra
using TightBindingLattice

@testset "basic" begin
    @testset "parser" begin
        using TightBindingLattice: parse_expr
        @test 1 == parse_expr(1)
        @test 1 == parse_expr("1")
        @test 1.5 == parse_expr("1.5")
        @test 1.5im == parse_expr("1.5im")
        @test 1.5im == parse_expr("1.5i")
        @test [1,2,3] == parse_expr("[1,2,3]")
        @test [1,2,3im] == parse_expr([1, "2", "3i"])
    end

    @testset "cleanup" begin
        using TightBindingLattice: cleanup_number
        @test cleanup_number(42, 1E-8) == 42
        @test cleanup_number(1.5 + 1E-12, 1E-8) == 1.5
        @test cleanup_number([1.0 + 1E-12, 2.0 - 1E-12, 3.0], 1E-8) == [1.0, 2.0, 3.0]
        @test cleanup_number(0.1234567, 1E-8) == 0.1234567
    end

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

            @test ExactLinearAlgebra.determinant(float.(mat)) == det(mat)
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

            @test maximum(abs.(inv(mat) - ExactLinearAlgebra.inverse(float.(mat)))) < sqrt(eps(Float64))
            @test maximum(abs.(inv(mat) - ExactLinearAlgebra.inverse(mat))) < sqrt(eps(Float64))
            for n in [5, 6, 7, 8]
                a = rand(-3:3, (n, n)) + 30*I
                inv1 = inv(a)
                inv2 = ExactLinearAlgebra.inverse(a)
                @test maximum(abs.(inv1 - inv2)) < sqrt(eps(Float64))
            end

        end
    end # testset ExactLinearAlgebra


end
