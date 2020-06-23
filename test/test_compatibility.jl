using Test
using TightBindingLattice


@testset "iscompatible" begin
    # hypercube and operation
    @testset "orthocube-operation" begin
        #hypercube = HypercubicLattice([2 0; 0 2])
        orthocube = OrthoCube([2 0; 0 2])
        @test_throws DimensionMismatch iscompatible(orthocube, TranslationOperation([1, 0, 0]))
        @test iscompatible(orthocube, TranslationOperation([1, 0]))
        
        @test_throws DimensionMismatch iscompatible(orthocube, PointOperation([1 0 0; 0 1 0; 0 0 0]))
        @test iscompatible(orthocube, PointOperation([0 1; 1 0]))
        @test !iscompatible(orthocube, PointOperation([0 2; 1 0]))
        @test !iscompatible(OrthoCube([3 0; 0 2]), PointOperation([0 1; 1 0]))
    end

    # hypercube and symmetry
    @testset "hypercube-symmetry" begin
        orthocube = OrthoCube([2 0; 0 2])
        @test  iscompatible(orthocube, TranslationSymmetry([2 0; 0 2]))
        @test  iscompatible(orthocube, TranslationSymmetry([2 2; 0 2]))
        @test !iscompatible(orthocube, TranslationSymmetry([1 0; 0 1]))
        @test !iscompatible(orthocube, TranslationSymmetry([2 0; 0 1]))

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(orthocube, psym)
        @test  iscompatible(OrthoCube([2 2; 0 2]), psym)
        @test !iscompatible(OrthoCube([3 0; 0 2]), psym)
    end
    
    # translation symmetry and point operation
    @testset "translation-point" begin
        @test  iscompatible(TranslationSymmetry([2 0; 0 2]), PointOperation([0 1; 1 0]))
        @test  iscompatible(TranslationSymmetry([2 2; 0 2]), PointOperation([0 1; 1 0]))
        @test !iscompatible(TranslationSymmetry([3 0; 0 2]), PointOperation([0 1; 1 0]))

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(TranslationSymmetry([2 0; 0 2]), psym)
        @test  iscompatible(TranslationSymmetry([2 2; 0 2]), psym)
        @test !iscompatible(TranslationSymmetry([3 0; 0 2]), psym)

        @test  iscompatible(psym, TranslationSymmetry([2 0; 0 2]))
        @test  iscompatible(psym, TranslationSymmetry([2 2; 0 2]))
        @test !iscompatible(psym, TranslationSymmetry([3 0; 0 2]))
    end


    @testset "translation-irrep-point-operation" begin
        # shape matching, momentum check
        tsym = TranslationSymmetry([3 0; 0 3])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test  iscompatible(tsym, irrep_index, PointOperation([0 1; 1  0])) == (k[1] == k[2])
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end
        tsym = TranslationSymmetry([3 3; 0 3])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test  iscompatible(tsym, irrep_index, PointOperation([0 1; 1  0])) == (k[1] == k[2])
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end

        # shape mismatch
        tsym = TranslationSymmetry([3 0; 0 2])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test !iscompatible(tsym, irrep_index, PointOperation([0 1; 1 0])) # not compatible with tsym
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(TranslationSymmetry([2 0; 0 2]), psym)
        @test  iscompatible(TranslationSymmetry([2 2; 0 2]), psym)
        @test !iscompatible(TranslationSymmetry([3 0; 0 2]), psym)
    end

    @testset "translation-irrep-point-symmetry" begin
        tsym = TranslationSymmetry([3 0; 0 3])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        tsym = TranslationSymmetry([4 0; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        tsym = TranslationSymmetry([4 4; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        # shape mismatch
        tsym = TranslationSymmetry([4 2; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test !iscompatible(tsym, irrep_index, psym)
            # 4mm already not compatible with lattice shape
        end
    end

    @testset "lattice-symmetry" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addsite!(unitcell, "X", FractCoord([0, 0], [0.5, 0.0]))
        
        lattice = make_lattice(unitcell, [2 0; 0 2])
        @test_throws DimensionMismatch iscompatible(lattice, TranslationOperation([1,2,3]))
        for i in -5:5, j in -5:5
            @test iscompatible(lattice, TranslationOperation([i, j]))
        end
        @test iscompatible(lattice, PointOperation([1 0; 0 -1]))
        @test !iscompatible(lattice, PointOperation([0 1; 1 0]))

        @test iscompatible(lattice, TranslationSymmetry([2 0; 0 2]))
        @test iscompatible(lattice, TranslationSymmetry([2 2; 0 2]))
        # ^- This is a choice. Whether to require equality, or just compatibility

        psym = project(PointSymmetryDatabase.get(7), [1 0 0; 0 1 0])
        @test iscompatible(lattice, psym)

        psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
        @test !iscompatible(lattice, psym)
    end
end

