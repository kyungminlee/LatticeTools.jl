using Test
using TightBindingLattice
using LinearAlgebra

@testset "embedding" begin
    include("Kagome.jl")

    kagome = make_kagome_lattice([4 0; 0 3])
    tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
    @test_throws ArgumentError embed(kagome.lattice, kagome.point_symmetry)
    @test_throws ArgumentError embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)
    
    kagome = make_kagome_lattice([2 0; 0 2])
    tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
    psym_embed = embed(kagome.lattice, kagome.point_symmetry)
    ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

    kagome = make_kagome_lattice([4 -2; 2 2])
    tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
    psym_embed = embed(kagome.lattice, kagome.point_symmetry)
    ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

    
    
end

