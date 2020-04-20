using Test
using TightBindingLattice

using LinearAlgebra


@testset "basic" begin
    @testset "round" begin
        for Ti in [Int, Int8, Int16, Int32, Int64]
            @test round(Ti,  1//2, RoundDown) == 0
            @test round(Ti, -1//2, RoundDown) == -1
            @test round(Ti,  1//2, RoundUp) == 1
            @test round(Ti, -1//2, RoundUp) == 0
            @test round(Ti,  1//2, RoundToZero) == 0
            @test round(Ti, -1//2, RoundToZero) == 0

            @test round(Ti,  1//1, RoundDown) == 1
            @test round(Ti, -1//1, RoundDown) == -1
            @test round(Ti,  1//1, RoundUp) == 1
            @test round(Ti, -1//1, RoundUp) == -1
            @test round(Ti,  1//1, RoundToZero) == 1
            @test round(Ti, -1//1, RoundToZero) == -1

            inf = 1 // 0
            @test_throws DivideError round(Ti,  inf, RoundDown)
            @test_throws DivideError round(Ti, -inf, RoundDown)
            @test_throws DivideError round(Ti,  inf, RoundUp)
            @test_throws DivideError round(Ti, -inf, RoundUp)
            @test_throws DivideError round(Ti,  inf, RoundToZero)
            @test_throws DivideError round(Ti, -inf, RoundToZero)
        end

        for Tf in [Float16, Float32, Float64]
            @test round(Tf,  1//2, RoundDown) == 0.0
            @test round(Tf, -1//2, RoundDown) == -1.0
            @test round(Tf,  1//2, RoundUp) == 1.0
            @test round(Tf, -1//2, RoundUp) == 0.0
            @test round(Tf,  1//2, RoundToZero) == 0.0
            @test round(Tf, -1//2, RoundToZero) == 0.0

            @test round(Tf,  1//1, RoundDown) == 1.0
            @test round(Tf, -1//1, RoundDown) == -1.0
            @test round(Tf,  1//1, RoundUp) == 1.0
            @test round(Tf, -1//1, RoundUp) == -1.0
            @test round(Tf,  1//1, RoundToZero) == 1.0
            @test round(Tf, -1//1, RoundToZero) == -1.0

            inf = 1 // 0
            @test round(Tf,  inf, RoundDown) == Inf
            @test round(Tf, -inf, RoundDown) == -Inf
            @test round(Tf,  inf, RoundUp) == Inf
            @test round(Tf, -inf, RoundUp) == -Inf
            @test round(Tf,  inf, RoundToZero) == Inf
            @test round(Tf, -inf, RoundToZero) == -Inf
        end
    end

    @testset "gcd" begin
        using TightBindingLattice: extended_gcd
        let x = 3, y = 0
            r, (a,b) = extended_gcd(x,y)
            @test r == 3
            @test a*x + b*y == r
        end
        let x = 0, y = -3
            r, (a,b) = extended_gcd(x,y)
            @test r == 3
            @test a*x + b*y == r
        end

        for sx in [-1,1], sy in [-1,1]
            let x = 3*sx, y = 4*sy
                r, (a,b) = extended_gcd(x,y)
                @test r == 1
                @test a*x + b*y == r
            end
            let x = 6*sx, y = 9*sy
                r, (a,b) = extended_gcd(x,y)
                @test r == 3
                @test a*x + b*y == r
            end
        end
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
