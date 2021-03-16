using Test
using LatticeTools

using LinearAlgebra


@testset "embedding" begin
    @testset "two-band" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
        addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
        addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
        @testset "square" begin
            lattice = make_lattice(unitcell, [2 0; 0 2])

            # |       |                 |       |
            # 6       8                 8       6
            # |       |         [1,0]   |       |
            # . - 5 - . - 7 -     =>    . - 7 - . - 5 -
            # |       |                 |       |
            # 2       4                 4       2
            # |       |                 |       |
            # o - 1 - . - 3 -           o - 3 - . - 1 -
            tsym = TranslationSymmetry(lattice)
            tsymbed = embed(lattice, tsym)
            translation_symmetry_embedding(lattice)


            let n = lowercase(symmetry_name(tsymbed))
                @test occursin("embed", n)
                @test occursin("translation", n)
                @test occursin("2", n) && occursin("0", n)
            end

            @test length(tsymbed) == 4
            @test length(collect(tsymbed)) == 4
            @test length(elements(tsymbed)) == 4
            @test eltype(tsymbed) == SitePermutation
            @test valtype(tsymbed) == SitePermutation
            @test eltype(typeof(tsymbed)) == SitePermutation
            @test valtype(typeof(tsymbed)) == SitePermutation
            @test elements(tsymbed)[1] == SitePermutation([1,2,3,4,5,6,7,8])
            @test elements(tsymbed)[2] == SitePermutation([3,4,1,2,7,8,5,6])
            @test element(tsymbed, 2) == SitePermutation([3,4,1,2,7,8,5,6])
            @test element(tsymbed, [1,2]) == [SitePermutation([1,2,3,4,5,6,7,8]), SitePermutation([3,4,1,2,7,8,5,6])]
            @test generator_indices(tsymbed) == [2, 3]
            @test generator_elements(tsymbed) == [element(tsymbed, 2) , element(tsymbed, 3)]

            # |       |                 |       |
            # 6       8                 3       7
            # |       |         C4      |       |
            # . - 5 - . - 7 -   =>      . - 8 - . - 4 -
            # |       |                 |       |
            # 2       4                 1       5
            # |       |                 |       |
            # o - 1 - . - 3 -           o - 6 - . - 2 -
            psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
            @test psym.hermann_mauguin == "4mm"
            psymbed = embed(lattice, psym)

            let n = lowercase(symmetry_name(psymbed))
                @test occursin("embed", n)
                @test occursin("point", n)
                @test occursin("4mm", n) || occursin("C<sub>4v</sub>", n)
            end

            @test eltype(psymbed) == SitePermutation
            @test valtype(psymbed) == SitePermutation

            idx_C4 = 3
            @test element_name(psym, idx_C4) == "4<sup>+</sup><sub>001</sub>"
            @test length(psymbed) == 8
            @test length(elements(psymbed)) == 8
            @test length(collect(psymbed)) == 8
            @test element(psymbed, idx_C4) == SitePermutation([2,3,6,7,4,1,8,5])

            @test iscompatible(tsymbed, psymbed)

            @test little_symmetry(tsymbed, psymbed).symmetry.hermann_mauguin == "4mm"
            @test little_symmetry(tsymbed, 1, psymbed).symmetry.hermann_mauguin == "4mm"
            @test little_symmetry(tsymbed, 2, psymbed).symmetry.hermann_mauguin == "mm2"
            @test little_symmetry(tsymbed, 3, psymbed).symmetry.hermann_mauguin == "mm2"
            @test little_symmetry(tsymbed, 4, psymbed).symmetry.hermann_mauguin == "4mm"

            @test little_group_elements(tsymbed, psymbed) == 1:group_order(psymbed)
            @test little_group_elements(tsymbed, 1, psymbed) == 1:group_order(psymbed)
            @test length(little_group_elements(tsymbed, 2, psymbed)) < group_order(psymbed)

            @test iscompatible(tsymbed, 1, psymbed)
            @test !iscompatible(tsymbed, 2, psymbed)
            @test !iscompatible(tsymbed, 3, psymbed)
            @test iscompatible(tsymbed, 4, psymbed)

            @test_throws ArgumentError embed(lattice, PointOperation([1 -1; 1 1])) # lattice not invariant under operation
            embed(make_lattice(unitcell, [1 0; 0 1]), psym) # lattice too small for faithful psym. Only a warning

            let lattice_large = make_lattice(unitcell, [4 0; 0 4])
                tsymbed_large = translation_symmetry_embedding(lattice_large)
                @test !iscompatible(tsymbed_large, psymbed)
                for irrep_index in 1:num_irreps(tsymbed_large)
                    @test !iscompatible(tsymbed_large, irrep_index, psymbed)
                end
                @test_throws ArgumentError little_group_elements(tsymbed_large, psymbed)
                @test_throws ArgumentError little_group_elements(tsymbed_large, 1, psymbed)
                @test_throws ArgumentError little_symmetry(tsymbed_large, psymbed)
                @test_throws ArgumentError little_symmetry(tsymbed_large, 1, psymbed)
            end
        end

        @testset "rectangular" begin
            lattice = make_lattice(unitcell, [3 0; 0 2])
            # tsym = TranslationSymmetry(lattice)
            # tsymbed = embed(lattice, tsym)
            psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
            @test_throws ArgumentError embed(lattice, psym)
            # @show little_symmetry(tsymbed, psymbed)
        end
    end

    @testset "kagome" begin
        include("Kagome.jl")
        # TODO: add tests
        kagome = make_kagome_lattice([4 0; 0 3])
        let
            kagome2 = make_kagome_lattice([4 0; 0 2])
            @test_throws ArgumentError embed(kagome.lattice, kagome2.translation_symmetry)
        end
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        @test_throws ArgumentError embed(kagome.lattice, kagome.point_symmetry)
        @test_throws ArgumentError embed(kagome.lattice, kagome.translation_symmetry ⋊ kagome.point_symmetry)
        @test_throws ArgumentError embed(kagome.lattice, kagome.translation_symmetry ⋊ˢ kagome.point_symmetry)

        # lattice too small for faithful embedding
        kagome = make_kagome_lattice([1 0; 0 1])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psymbed = embed(kagome.lattice, kagome.point_symmetry)
        ssymbed = embed(kagome.lattice, kagome.translation_symmetry ⋊ kagome.point_symmetry)
        ssymbed2 = embed(kagome.lattice, kagome.translation_symmetry ⋊ˢ kagome.point_symmetry)

        kagome = make_kagome_lattice([2 0; 0 2])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psym_embed = embed(kagome.lattice, kagome.point_symmetry)
        ssym_embed = embed(kagome.lattice, kagome.translation_symmetry ⋊ kagome.point_symmetry)
        ssym_embed2 = embed(kagome.lattice, kagome.translation_symmetry ⋊ˢ kagome.point_symmetry)

        @test group_order(tsym_embed) == group_order(kagome.translation_symmetry)
        @test group_order(psym_embed) == group_order(kagome.point_symmetry)
        @test group_order(ssym_embed) == group_order(kagome.translation_symmetry) * group_order(kagome.point_symmetry)
        @test group_order(ssym_embed2) == group_order(kagome.translation_symmetry) * group_order(kagome.point_symmetry)

        kagome = make_kagome_lattice([4 -2; 2 2])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psym_embed = embed(kagome.lattice, kagome.point_symmetry)
        ssym_embed = embed(kagome.lattice, kagome.translation_symmetry ⋊ kagome.point_symmetry)
    end

    # @testset "kagome-strong" begin
    #     include("Kagome.jl")
    #     kagome = make_kagome_lattice([4 -2; 2 2])
    #     tsymbed = embed(kagome.lattice, kagome.translation_symmetry)
    #     psymbed = little_symmetry_strong(tsymbed, embed(kagome.lattice, kagome.point_symmetry))
    #     @show symmetry_name(tsymbed)
    #     @show symmetry_name(psymbed)
    # end
end
