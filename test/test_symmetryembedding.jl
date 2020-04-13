using Test
using TightBindingLattice
using LinearAlgebra

@testset "embedding" begin
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
