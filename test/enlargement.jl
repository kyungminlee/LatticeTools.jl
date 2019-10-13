using Test
using LinearAlgebra
using TightBindingLattice

@testset "enlargement.jl" begin
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
      @test_throws ArgumentError ExactLinearAlgebra.inverse(ones(Int, 0,0))
      @test ExactLinearAlgebra.inverse(7*ones(Int, 1, 1)) == ones(Int, 1, 1)//7

      @test maximum(abs.(inv(mat) - ExactLinearAlgebra.inverse(mat))) < sqrt(eps(Float64))
      for n in [5, 6, 7, 8]
        a = rand(-3:3, (n, n))
        inv1 = inv(a)
        inv2 = ExactLinearAlgebra.inverse(a)
        @test maximum(abs.(inv1 - inv2)) < sqrt(eps(Float64))
      end

    end
  end # testset ExactLinearAlgebra



end # testset enlargement
