#=
using Test
using LatticeTools


@testset "Permutation" begin
    @test_throws ArgumentError Permutation([1,2,4])
    @test_throws OverflowError Permutation([mod(x, 4096)+1 for x in 1:4096])
    p0 = Permutation([1,2,3,4])
    p1 = Permutation([2,3,4,1])
    p2 = Permutation([3,4,1,2])
    p3 = Permutation([4,1,2,3])

    @test p1 * p2 == p3
    @test p1 != p3
    @test p1^0 == p0
    @test p1^1 == p1
    @test p1^2 == p2
    @test p1^3 == p3

    @test p0.order == 1
    @test p1.order == 4
    @test p2.order == 2
    @test p3.order == 4

    @test inv(p0) == p0
    @test inv(p1) == p3
    @test inv(p2) == p2
    @test inv(p3) == p1

    @test p0 == Permutation([1,2,3,4])
    @test isless(p0, p1)
    @test isless(p1, p3)

    @test p0(1) == 1 && p0(2) == 2 && p0(3) == 3 && p0(4) == 4
    @test p1(1) == 2 && p1(2) == 3 && p1(3) == 4 && p1(4) == 1

    @test_throws ArgumentError Permutation([1,2,3,4]) * Permutation([1,2,3,4,5])
    @test hash(Permutation(Int[1,2,3,4])) != hash(Int[1,2,3,4])

    @test generate_group(p1) == Set([p0, p1, p2, p3])
    @test generate_group(p2) == Set([p0, p2])
    @test generate_group(p1, p2) == Set([p0, p1, p2, p3])
end
=#