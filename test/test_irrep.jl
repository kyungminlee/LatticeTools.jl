using Test
using LatticeTools


@testset "irrep" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
    addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
    addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

    lattice = make_lattice(unitcell, [4 0; 0 4])
    tsym = TranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

    tsymbed = embed(lattice, tsym)
    psymbed = embed(lattice, psym)

    @testset "translation" begin
        @test_throws ArgumentError IrrepComponent(tsym, 99999)
        @test_throws ArgumentError IrrepComponent(tsym, 1, 2)
        @test IrrepComponent(tsym, 1, 1) == IrrepComponent(tsym, 1, 1)
        @test IrrepComponent(tsym, 1, 1) != IrrepComponent(tsym, 2, 1)

        for tsym_irrep in 1:num_irreps(tsym)
            tsic = IrrepComponent(tsym, tsym_irrep)
            tsicbed = IrrepComponent(tsymbed, tsym_irrep)

            @test group_order(tsic) == group_order(tsym)
            psym_little = little_symmetry(tsic, psym)
            psymbed_little = little_symmetry(tsicbed, psymbed)
            @test psym_little.group == symmetry(psymbed_little).group

            k = tsym.coordinates[tsym_irrep]
            @test iscompatible(tsym, tsym_irrep, psym) == (k in [[0,0], [2,2]])
            @test iscompatible(tsym, tsym_irrep, psym_little)

            @test iscompatible(tsic, psym) == (k in [[0,0], [2,2]])
            @test iscompatible(tsic, psym_little)

            @test iscompatible(tsicbed, psymbed) == (k in [[0,0], [2,2]])
            @test iscompatible(tsicbed, psymbed_little)

            @test little_group(tsic, psym) == little_group(tsicbed, psymbed)

            lge1 = little_group_elements(tsic, psym)
            lge2 = little_group_elements(tsicbed, psymbed)
            @test lge1 == lge2
            lg_matrep = psym.matrix_representations[lge1]
            @test !isnothing(group_isomorphism(little_group(tsym, tsym_irrep, psym),
                                  FiniteGroup(group_multiplication_table(lg_matrep))))
            let psic = IrrepComponent(psym, 1, 1)
                if k in [[0,0], [2,2]]
                    SymmorphicIrrepComponent(tsic, psic)
                else
                    @test_throws ArgumentError SymmorphicIrrepComponent(tsic, psic)
                end
            end
        end # for tsym_irrep

        for tsic in get_irrep_components(tsym)
            k = tsym.coordinates[tsic.irrep_index]
            psym_little = little_symmetry(tsic, psym)

            @test iscompatible(tsic, psym) == (k in [[0,0], [2,2]])
            @test iscompatible(tsic, psym_little)
            lg_matrep = psym.matrix_representations[little_group_elements(tsic, psym)]
            @test !isnothing(group_isomorphism(little_group(tsic, psym),
                                FiniteGroup(group_multiplication_table(lg_matrep))))
            let psic = IrrepComponent(psym, 1, 1)
                if k in [[0,0], [2,2]]
                    SymmorphicIrrepComponent(tsic, psic)
                else
                    @test_throws ArgumentError SymmorphicIrrepComponent(tsic, psic)
                end
            end
        end
    end

    @testset "point" begin
        @test_throws ArgumentError IrrepComponent(psym, 99999, 1)
        @test_throws ArgumentError IrrepComponent(psym, 1, 10)
        @test length(collect(get_irrep_components(psym))) == 1 + 1 + 1 + 1 + 2
        @test IrrepComponent(psym, 1, 1) == IrrepComponent(psym, 1, 1)
        @test IrrepComponent(psym, 1, 1) != IrrepComponent(psym, 2, 1)

        for psic in get_irrep_components(psym)
            @test group_order(psic) == group_order(psym)

            permphase = collect(get_irrep_iterator(psic))
            @test length(permphase) == group_order(psym)

            perms = [p for (p, ϕ) in permphase]
            @test ishomomorphic(psym.group, perms)

            mtab = group_multiplication_table(perms)
            @test mtab == psym.group.multiplication_table
        end
    end

    @testset "space" begin
        ssics1 = []
        for tsic in get_irrep_components(tsym)
            psym_little = little_symmetry(tsic, psym)
            for psic in get_irrep_components(psym_little)
                push!(ssics1, (tsic, psic))
            end
        end
        ssics2 = []
        for ssic in get_irrep_components(tsym ⋊ psym)
            push!(ssics2, ssic)
        end
        @test length(ssics1) == length(ssics2)
        for ((tsic, psic), ssic) in zip(ssics1, ssics2)
            @test group_order(ssic) == group_order(tsic) * group_order(psic)

            @test tsic.symmetry.hypercube == ssic.normal.symmetry.hypercube
            @test psic.symmetry.hermann_mauguin == ssic.rest.symmetry.hermann_mauguin
            @test tsic.irrep_index == ssic.normal.irrep_index
            @test psic.irrep_index == ssic.rest.irrep_index
            @test tsic.irrep_component == ssic.normal.irrep_component
            @test psic.irrep_component == ssic.rest.irrep_component

            @test sum(1 for x in collect(get_irrep_iterator(ssic))) == group_order(tsic) * group_order(psic)
        end
    end

end # testset little_symmetry
