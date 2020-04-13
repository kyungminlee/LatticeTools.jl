using Test
using TightBindingLattice
using LinearAlgebra

@testset "embedding" begin


    @testset "two-band" begin

        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
        addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
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
        # @show tsymbed.elements[2]
        @test length(elements(tsymbed)) == 4
        @test element(tsymbed, 2) == SitePermutation([3,4,1,2,7,8,5,6])

        # |       |                 |       |
        # 6       8                 3       7
        # |       |         C4      |       |
        # . - 5 - . - 7 -   =>      . - 8 - . - 4 -
        # |       |                 |       |
        # 2       4                 1       5
        # |       |                 |       |
        # o - 1 - . - 3 -           o - 6 - . - 2 -
        psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
        idx_C4 = 3
        @test element_name(psym, idx_C4) == "4<sup>+</sup><sub>001</sub>"
        psymbed = embed(lattice, psym)
        @test length(elements(psymbed)) == 8
        @test element(psymbed, idx_C4) == SitePermutation([2,3,6,7,4,1,8,5])
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
        @test_throws ArgumentError embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

        kagome = make_kagome_lattice([1 0; 0 1])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psym_embed = embed(kagome.lattice, kagome.point_symmetry)
        ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

        kagome = make_kagome_lattice([2 0; 0 2])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psym_embed = embed(kagome.lattice, kagome.point_symmetry)
        ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

        @test group_order(tsym_embed) == group_order(kagome.translation_symmetry)
        @test group_order(psym_embed) == group_order(kagome.point_symmetry)
        #@test group_order(ssym_embed) == group_order(kagome.translation_symmetry) * group_order(kagome.point_symmetry)


        kagome = make_kagome_lattice([4 -2; 2 2])
        tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
        psym_embed = embed(kagome.lattice, kagome.point_symmetry)
        ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)
    end
end
