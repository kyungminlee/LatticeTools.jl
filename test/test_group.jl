using Test

using LinearAlgebra
using YAML

@testset "Group" begin
    @testset "FiniteGroup-Abelian" begin
        @test_throws ArgumentError FiniteGroup([1 1 1; 1 1 1])
        @test_throws ArgumentError FiniteGroup([1 1 1; 1 1 1; 1 1 1])

        # Example: Z₃
        mtab = [1 2 3;
                2 3 1;
                3 1 2]
        group = FiniteGroup(mtab)

        @test group_order(group) == 3
        @test group_order(group, 1) == 1
        @test group_order(group, 2) == 3
        @test group_order(group, 3) == 3
        @test isabelian(group)
        @test group_multiplication_table(group) == mtab

        for i in 1:3, j in 1:3
            @test group_product(group, i, j) == mtab[i, j]
        end
        @test group_product(group, 2, BitSet([1,2])) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), 2) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), BitSet([1,2])) == BitSet([1,2,3])

        @test generate_subgroup(group, 1) == BitSet([1])
        @test generate_subgroup(group, 2) == BitSet([1,2,3])
        @test generate_subgroup(group, [1,2]) == BitSet([1,2,3])

        @test issubgroup(group, Set([1]))
        @test !issubgroup(group, Set([1,2]))
        @test minimal_generating_set(group) == [2]

        @test group_multiplication_table([[1 0; 0 1], [1 0; 0 -1]]) == [1 2; 2 1]
    end

    @testset "FiniteGroup-Nonabelian" begin
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


        ϕ = [1, 3, 4, 5, 2, 6] # group isomorphism
        mtab1 = group_multiplication_table(group)
        mtab2 = zeros(Int, (6,6))
        for x in 1:6, y in 1:6
            # ϕ(x)⋅ϕ(y) = ϕ(x⋅y)
            mtab2[ϕ[x], ϕ[y]] = ϕ[mtab1[x,y]]
        end
        group2 = FiniteGroup(mtab2)

        # ϕ: group  ->  group2
        #       x  |->  ϕ(x)
        ϕ2 = group_isomorphism(group, group2)
        mtab3 = zeros(Int, (6,6))
        for x in 1:6, y in 1:6
            # ϕ(x)⋅ϕ(y) = ϕ(x⋅y)
            mtab3[ϕ2[x], ϕ2[y]] = ϕ2[mtab1[x,y]]
        end
        @test !isnothing(group_isomorphism(group2, FiniteGroup(mtab3)))  # ϕ and ϕ2 are equivalent

        #@show mtab2
        #@show group_isomorphism(group, group2) # finds ϕ
        #@show group_isomorphism(group2, group) # finds ϕ⁻¹

    end

    @testset "group isomorphism" begin
        # 4
        group1 = FiniteGroup([1 2 3 4;
                              2 1 4 3;
                              3 4 2 1;
                              4 3 1 2])
        group1p= FiniteGroup([1 2 3 4;
                              2 3 4 1;
                              3 4 1 2;
                              4 1 2 3])
        @test !isnothing(group_isomorphism(group1, group1p))
        # 2/m
        group2 = FiniteGroup([1 2 3 4;
                              2 1 4 3;
                              3 4 1 2;
                              4 3 2 1])
        @test isnothing(group_isomorphism(group1, group2))
    end # @testset "group isomorphism" begin
end
