using Test
using LatticeTools

@testset "linpath" begin
    @test_throws DimensionMismatch linpath([0.0, 1.0], [0.0, 1.0, 2.0])
    @test_throws ArgumentError linpath([0.0, 1.0])

    r = linpath([0.0, 0.0], [1.0, 0.0], [1.0, 1.0])
    @test size(r) == (2, 100*2+1)
    @test r[1,:] ≈ vcat(0:0.01:0.99, [1.0 for i in 1:100], [1.0])
    @test r[2,:] ≈ vcat([0.0 for i in 1:100], 0:0.01:0.99, [1.0])

    r = linpath([0.0, 0.0], [1.0, 0.0], [1.0, 1.0]; nseg=10)
    @test size(r) == (2, 10*2+1)
    @test r[1,:] ≈ vcat(0:0.1:0.9, [1.0 for i in 1:10], [1.0])
    @test r[2,:] ≈ vcat([0.0 for i in 1:10], 0:0.1:0.9, [1.0])
end
