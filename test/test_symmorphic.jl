using Test
using TightBindingLattice

@testset "symmorphic" begin
    tsym = TranslationSymmetry([4 0; 0 4])
    psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
    @test_throws ArgumentError TranslationSymmetry([4 0; 0 3]) ⋊ psym
    @test_throws ArgumentError psym ⋊ tsym
    ssym = tsym ⋊ psym

    @test eltype(ssym) == SpaceOperation{Int, Int}
    @test valtype(ssym) == SpaceOperation{Int, Int}
    @test element(ssym, 2) == TranslationOperation([1, 0])  # equality between different types
    let els = elements(ssym)
        @test length(els) == 4 * 4 * group_order(psym)
        @test typeof(first(els)) == SpaceOperation{Int, Int}
    end

    @test element_name(ssym, 2) == "{ [1, 0] | 1 }"
    let els = element_names(ssym)
        @test length(els) == 4 * 4 * group_order(psym)
        @test typeof(first(els)) <: AbstractString
    end

    els = elements(ssym)
    let product = symmetry_product(ssym)
        @test all(let
                      z1 = product(x, y)
                      z2 = x * y
                      z3 = SpaceOperation(z2.matrix, tsym.orthocube.wrap(z2.displacement)[2])
                      z1 == z3
                  end for x in els, y in els)
    end
    element_lookup = Dict(x => i for (i, x) in enumerate(elements(ssym)))
    
    @test group_order(ssym) == group_order(psym) * group_order(tsym)

    let
        mtab = zeros(Int, group_order(ssym), group_order(ssym))
        let product = symmetry_product(ssym)
            for (i, xi) in enumerate(els), (j, xj) in enumerate(els)
                xk = product(xi, xj)
                k = element_lookup[xk]
                mtab[i,j] = k
            end
        end
        G = group(ssym)
        @test G.multiplication_table == mtab
        @test group_multiplication_table(ssym) == mtab
    end

    gen_inds = generator_indices(ssym)
    gen_els = generator_elements(ssym)
    @test length(gen_inds) == length(generator_indices(tsym)) + length(generator_indices(psym))
    @test length(gen_els) == length(generator_indices(tsym)) + length(generator_indices(psym))

    # get_irrep_components(ssym)

    @testset "embedding" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(unitcell, "A", FractCoord([0, 0], [0.5, 0.0]))
        addorbital!(unitcell, "B", FractCoord([0, 0], [0.0, 0.5]))
        lattice = make_lattice(unitcell, [4 0; 0 4])
        tsymbed = embed(lattice, tsym)
        psymbed = embed(lattice, psym)

        ssym = tsym ⋊ psym

        ssymbed1 = embed(lattice, ssym)
        ssymbed2 = tsymbed ⋊ psymbed
        ssymbed3 = SymmorphicSymmetryEmbedding(tsymbed, psymbed)
        @test_throws ArgumentError SymmorphicSymmetryEmbedding(psymbed, tsymbed)
        @test_throws ArgumentError SymmorphicSymmetryEmbedding(make_lattice(unitcell, [4 0; 0 3]), ssym)
        
        let lattice2 = make_lattice(unitcell, [4 0; 0 3])
            tsymbed2 = translation_symmetry_embedding(lattice2)
            @test_throws ArgumentError tsymbed2 ⋊ psymbed
        end

        @test typeof(ssymbed1) == typeof(ssymbed2)
        @test eltype(ssymbed1) == SitePermutation
        @test eltype(ssymbed2) == SitePermutation
        @test valtype(ssymbed1) == SitePermutation
        @test valtype(ssymbed2) == SitePermutation

        @test ssymbed1 == ssymbed2
        @test ssymbed1 == ssymbed3
        
        @test ssymbed1.lattice == lattice
        @test ssymbed1.normal.symmetry == tsymbed.symmetry
        @test ssymbed1.rest.symmetry == psymbed.symmetry

        @test ssymbed2.lattice == lattice
        @test ssymbed2.normal.symmetry == tsymbed.symmetry
        @test ssymbed2.rest.symmetry == psymbed.symmetry

        tsic1 = first(get_irrep_components(tsymbed))
        psic1 = first(get_irrep_components(psymbed))
        SymmorphicIrrepComponent(tsic1, psic1)

        let lattice2 = make_lattice(unitcell, [4 0; 0 3])
            tsymbed2 = translation_symmetry_embedding(lattice2)
            tsic2 = first(get_irrep_components(tsymbed2))
            @test_throws ArgumentError SymmorphicIrrepComponent(tsic2, psic1)
        end

    end
end