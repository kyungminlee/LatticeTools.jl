using Test
using TightBindingLattice

using Combinatorics
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
        @test group_order(group, [1,2,3]) == [1,3,3]
        @test period_length(group, 1) == 1
        @test period_length(group, 2) == 3
        @test period_length(group, 3) == 3
        @test period_length(group, [1,2,3]) == [1,3,3]
        @test isabelian(group)
        @test group_multiplication_table(group) == mtab
        @test element(group, 2) == 2
        @test element(group, 1:2) == 1:2
        @test elements(group) == 1:3
        @test element_name(group, 2) == "2"
        @test element_name(group, 1:2) == ["1", "2"]
        @test element_names(group) == ["1", "2", "3"]
        @test_throws BoundsError element(group, 5)
        @test_throws BoundsError element_name(group, 5)

        gp = group_product(group)
        for i in 1:3, j in 1:3
            @test group_product(group, i, j) == mtab[i, j]
            @test gp(i, j) == mtab[i, j]
        end
        @test group_product(group, 2, BitSet([1,2])) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), 2) == BitSet([2,3])
        @test group_product(group, BitSet([1,2]), BitSet([1,2])) == BitSet([1,2,3])

        @test gp(2, BitSet([1,2])) == BitSet([2,3])
        @test gp(BitSet([1,2]), 2) == BitSet([2,3])
        @test gp(BitSet([1,2]), BitSet([1,2])) == BitSet([1,2,3])


        @test group_inverse(group, 1) == 1
        @test group_inverse(group, 2) == 3
        @test group_inverse(group, 3) == 2
        @test group_inverse(group, [1,2]) == [1,3]

        ginv = group_inverse(group)
        @test ginv(1) == 1
        @test ginv(2) == 3
        @test ginv(3) == 2
        @test ginv([1,2]) == [1,3]

        @test generate_subgroup(group, 1) == BitSet([1])
        @test generate_subgroup(group, 2) == BitSet([1,2,3])
        @test generate_subgroup(group, [1,2]) == BitSet([1,2,3])

        @test issubgroup(group, Set([1]))
        @test !issubgroup(group, Set([1,2]))
        @test minimal_generating_set(group) == [2]

        @test group_multiplication_table([[1 0; 0 1], [1 0; 0 -1]]) == [1 2; 2 1]

        @test ishomomorphic(group, 1:3; product=gp)
        @test !ishomomorphic(group, 1:2; product=gp)
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

        # ϕ: group  →  group2
        #       x   ↦  ϕ(x)
        ϕ2 = group_isomorphism(group, group2)
        mtab3 = zeros(Int, (6,6))
        for x in 1:6, y in 1:6
            # ϕ(x)⋅ϕ(y) = ϕ(x⋅y)
            mtab3[ϕ2[x], ϕ2[y]] = ϕ2[mtab1[x,y]]
        end
        @test !isnothing(group_isomorphism(group2, FiniteGroup(mtab3)))  # ϕ and ϕ2 are equivalent

        for ϕ in permutations(2:6)
            ϕ = vcat([1], ϕ)
            mtab1 = group_multiplication_table(group)
            mtab2 = zeros(Int, (6,6))
            for x in 1:6, y in 1:6
                # ϕ(x)⋅ϕ(y) = ϕ(x⋅y)
                mtab2[ϕ[x], ϕ[y]] = ϕ[mtab1[x,y]]
            end
            group2 = FiniteGroup(mtab2)
            @test !isnothing(group_isomorphism(group, group2))
        end
    end

    @testset "minimal_generating_set" begin
        for i in 1:32
            psym = PointSymmetryDatabase.get(i)
            group = psym.group
            mgs = minimal_generating_set(group)
            @test length(mgs) <= length(psym.generators)
            @test generate_subgroup(group, mgs) == BitSet(1:group_order(group))
        end
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

    @testset "group isomorphism 2" begin
        # D2h
        group1 = FiniteGroup([1 2 3 4 5 6 7 8;
                              2 1 4 3 6 5 8 7;
                              3 4 1 2 7 8 5 6;
                              4 3 2 1 8 7 6 5;
                              5 6 7 8 1 2 3 4;
                              6 5 8 7 2 1 4 3;
                              7 8 5 6 3 4 1 2;
                              8 7 6 5 4 3 2 1])
        group2 = FiniteGroup([1 2 3 4 5 6 7 8;
                              2 1 4 3 6 5 8 7;
                              3 4 2 1 7 8 6 5;
                              4 3 1 2 8 7 5 6;
                              5 6 7 8 1 2 3 4;
                              6 5 8 7 2 1 4 3;
                              7 8 6 5 3 4 2 1;
                              8 7 5 6 4 3 1 2])
        @test isnothing(group_isomorphism(group1, group2))
    end

    @testset "group isomorphism 3" begin
        tsym1 = TranslationSymmetry([3 0; 0 3])
        tsym2 = TranslationSymmetry([3 1; 0 3])
        @test isnothing(group_isomorphism(tsym1.group, tsym2.group))
    end

    # TODO
    # Are there two non-isomorphic groups with
    # - Same conjugacy classes
    # - Same period lengths
end
