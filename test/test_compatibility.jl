using Test
using LatticeTools


@testset "iscompatible" begin
    # hypercube and operation
    @testset "hypercube-operation" begin
        #hypercube = HypercubicLattice([2 0; 0 2])
        hypercube = Hypercube([2 0; 0 2])
        @test_throws DimensionMismatch iscompatible(hypercube, TranslationOperation([1, 0, 0]))
        @test iscompatible(hypercube, TranslationOperation([1, 0]))

        @test_throws DimensionMismatch iscompatible(hypercube, PointOperation([1 0 0; 0 1 0; 0 0 0]))
        @test iscompatible(hypercube, PointOperation([0 1; 1 0]))
        @test !iscompatible(hypercube, PointOperation([0 2; 1 0]))
        @test !iscompatible(Hypercube([3 0; 0 2]), PointOperation([0 1; 1 0]))
    end

    # hypercube and symmetry
    @testset "hypercube-symmetry" begin
        hypercube = Hypercube([2 0; 0 2])
        @test  iscompatible(hypercube, FiniteTranslationSymmetry([2 0; 0 2]))
        @test  iscompatible(hypercube, FiniteTranslationSymmetry([2 2; 0 2]))
        @test !iscompatible(hypercube, FiniteTranslationSymmetry([1 0; 0 1]))
        @test !iscompatible(hypercube, FiniteTranslationSymmetry([2 0; 0 1]))

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(hypercube, psym)
        @test  iscompatible(Hypercube([2 2; 0 2]), psym)
        @test !iscompatible(Hypercube([3 0; 0 2]), psym)
    end

    # translation symmetry and point operation
    @testset "translation-point" begin
        @test  iscompatible(FiniteTranslationSymmetry([2 0; 0 2]), PointOperation([0 1; 1 0]))
        @test  iscompatible(FiniteTranslationSymmetry([2 2; 0 2]), PointOperation([0 1; 1 0]))
        @test !iscompatible(FiniteTranslationSymmetry([3 0; 0 2]), PointOperation([0 1; 1 0]))

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(FiniteTranslationSymmetry([2 0; 0 2]), psym)
        @test  iscompatible(FiniteTranslationSymmetry([2 2; 0 2]), psym)
        @test !iscompatible(FiniteTranslationSymmetry([3 0; 0 2]), psym)

        @test  iscompatible(psym, FiniteTranslationSymmetry([2 0; 0 2]))
        @test  iscompatible(psym, FiniteTranslationSymmetry([2 2; 0 2]))
        @test !iscompatible(psym, FiniteTranslationSymmetry([3 0; 0 2]))
    end


    @testset "translation-irrep-point-operation" begin
        # shape matching, momentum check
        tsym = FiniteTranslationSymmetry([3 0; 0 3])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test  iscompatible(tsym, irrep_index, PointOperation([0 1; 1  0])) == (k[1] == k[2])
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end
        tsym = FiniteTranslationSymmetry([3 3; 0 3])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test  iscompatible(tsym, irrep_index, PointOperation([0 1; 1  0])) == (k[1] == k[2])
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end

        # shape mismatch
        tsym = FiniteTranslationSymmetry([3 0; 0 2])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test !iscompatible(tsym, irrep_index, PointOperation([0 1; 1 0])) # not compatible with tsym
            @test  iscompatible(tsym, irrep_index, PointOperation([1 0; 0 -1])) == (k[2] == 0 || k[2] == 1//2)
        end

        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        @test  iscompatible(FiniteTranslationSymmetry([2 0; 0 2]), psym)
        @test  iscompatible(FiniteTranslationSymmetry([2 2; 0 2]), psym)
        @test !iscompatible(FiniteTranslationSymmetry([3 0; 0 2]), psym)
    end

    @testset "translation-irrep-point-symmetry" begin
        tsym = FiniteTranslationSymmetry([3 0; 0 3])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        tsym = FiniteTranslationSymmetry([4 0; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        tsym = FiniteTranslationSymmetry([4 4; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test iscompatible(tsym, irrep_index, psym) == (k in [[0//1, 0//1], [1//2, 1//2]])
            # 4mm only at Gamma and M points
        end

        # shape mismatch
        tsym = FiniteTranslationSymmetry([4 2; 0 4])
        psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            @test !iscompatible(tsym, irrep_index, psym)
            # 4mm already not compatible with lattice shape
        end
    end

    @testset "lattice-symmetry" begin
        unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
        addsite!(unitcell, "X", FractCoord([0, 0], [0.5, 0.0]))

        lattice = make_lattice(unitcell, [2 0; 0 2])
        @test_throws DimensionMismatch iscompatible(lattice, TranslationOperation([1,2,3]))
        for i in -5:5, j in -5:5
            @test iscompatible(lattice, TranslationOperation([i, j]))
        end
        @test iscompatible(lattice, PointOperation([1 0; 0 -1]))
        @test !iscompatible(lattice, PointOperation([0 1; 1 0]))

        @test iscompatible(lattice, FiniteTranslationSymmetry([2 0; 0 2]))
        @test iscompatible(lattice, FiniteTranslationSymmetry([2 2; 0 2]))
        # ^- This is a choice. Whether to require equality, or just compatibility

        psym = project(PointSymmetryDatabase.get(7), [1 0 0; 0 1 0])
        @test iscompatible(lattice, psym)

        psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
        @test !iscompatible(lattice, psym)
    end
end
