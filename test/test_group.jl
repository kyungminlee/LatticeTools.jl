using Test

using LinearAlgebra
using YAML

@testset "Group" begin
    @testset "FiniteAbelianGroup" begin
        @test_throws ArgumentError FiniteAbelianGroup([1 1 1; 1 1 1])
        @test_throws ArgumentError FiniteAbelianGroup([1 1 1; 1 1 1; 1 1 1])

        # Example: Z₃
        mtab = [1 2 3;
                2 3 1;
                3 1 2]
        group = FiniteAbelianGroup(mtab)

        @test order(group) == 3
        @test n_elements(group) == 3
        @test is_abelian(group)

        for i in 1:3, j in 1:3
            @test group_product(group, i, j) == mtab[i, j]
        end
        @test group_product(group, 2, BitSet([1,2])) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), 2) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), BitSet([1,2])) == BitSet([1,2,3])

        @test generate_subgroup(group, 1) == BitSet([1])
        @test generate_subgroup(group, 2) == BitSet([1,2,3])

        @test is_subgroup(group, Set([1]))
        @test !is_subgroup(group, Set([1,2]))
        @test minimal_generating_set(group) == [2]
    end

    @testset "FiniteGroup" begin
        # C3v, element#2 = C3, element#4 = σᵥₐ
        # (non-abelian)
        group = FiniteGroup([1 2 3 4 5 6;
                              2 3 1 6 4 5;
                              3 1 2 5 6 4;
                              4 5 6 1 2 3;
                              5 6 4 3 1 2;
                              6 4 5 2 3 1])
        generators = minimal_generating_set(group)
        @test generators == [2, 4]
        @test generate_subgroup(group, generators) == BitSet(1:6) # completely generates
    end
end
