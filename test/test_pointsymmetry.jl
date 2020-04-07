using Test
using LinearAlgebra
using YAML

using TightBindingLattice: simplify_name

@testset "PointSymmetry" begin
    @testset "simplify_name" begin
        @test simplify_name("1") == "1"
        @test simplify_name("-1") == "-1"
        @test simplify_name("2<sub>011</sub>") == "2"
        @test simplify_name("3<sup>+</sup><sub>011</sub>") == "3"
        @test simplify_name("m<sub>01</sub>") == "m"
    end

    @testset "constructor failure" begin
        group = FiniteGroup([1 2; 2 1])
        generators = [2]
        conjugacy_classes = [[1], [2]]
        character_table = [1 1; 1 -1]
        irreps = [
                [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                [ones(ComplexF64, 1, 1), -ones(ComplexF64, 1, 1)],
        ]
        element_names = ["1", "-1"]
        matrix_representations = [[1 0; 0 1], [-1 0; 0 -1]]
        hermann_mauguinn = "-1"
        schoenflies = "C<sub>i</sub>"
        PointSymmetry(group, generators,
                      conjugacy_classes, character_table, irreps,
                      element_names, matrix_representations, hermann_mauguinn, schoenflies)

        let generators = [1]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let conjugacy_classes = [[1], [1,2]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let character_table = [1 2 3; 4 5 6; 7 8 9]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [ones(ComplexF64, 1, 1), -ones(ComplexF64, 2, 2)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [2*ones(ComplexF64, 1, 1), -ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let element_names = ["1", "2", "3"]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let matrix_representations = [[1 0; 0 1], [-1 0 0; 0 -1 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let matrix_representations = [[1 0 0; 1 0 0], [-1 0 0; -1 0 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let matrix_representations = [[1 0 0; 1 0 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
        let matrix_representations = [[1 0 ; 0 1], [1 0; 0 1]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                                                                                              conjugacy_classes, character_table, irreps,
                                                                                              element_names, matrix_representations, hermann_mauguinn, schoenflies)
        end
    end

    # D_4
    file_path = abspath(@__DIR__, "..", "data", "PointGroup3D", "PointGroup3D-12.yaml")
    data_yaml = YAML.load_file(file_path)

    let data_yaml_2 = deepcopy(data_yaml)
        data_yaml_2["Generators"] = [1]
        @test_throws ArgumentError read_point_symmetry(data_yaml_2)
    end

    psym = read_point_symmetry(data_yaml)
    ord_group = group_order(psym.group)

    @testset "elements" begin
        @test group_order(psym) == 8
        for i in 1:ord_group
            @test element_names(psym)[i] == element_name(psym, i)
        end
    end

    @testset "generators" begin
        # Generators completely generate the group
        @test generate_subgroup(psym.group, psym.generators) == BitSet(1:ord_group)
        @test prod(psym.group.period_lengths[x] for x in psym.generators) != ord_group # nonabelian
    end

    @testset "multiplication table" begin
        matrep_lookup = Dict(m=>i for (i, m) in enumerate(psym.matrix_representations))
        matrep_mtab = zeros(Int, (ord_group, ord_group))
        for i1 in 1:ord_group, i2 in 1:ord_group
            m1 = psym.matrix_representations[i1]
            m2 = psym.matrix_representations[i2]
            m3 = m1 * m2
            i3 = matrep_lookup[m3]
            matrep_mtab[i1, i2] = i3
        end
        @test matrep_mtab == psym.group.multiplication_table
        @test matrep_mtab == group_multiplication_table(psym)
        @test matrep_mtab == group_multiplication_table(psym.group)
    end


    @testset "irreps" begin
        @test num_irreps(psym) == 5
        let χ = [1  1  1  1  1;
                          1  1  1 -1 -1;
                          1  1 -1 -1  1;
                          1  1 -1  1 -1;
                          2 -2  0  0  0]
            @test character_table(psym) ≈ χ
            @test size(character_table(psym)) == (5, 5)
        end

        @test length(irreps(psym)) == 5
        for i in 1:5
            @test irreps(psym)[i] == irrep(psym, i)
        end

        for (idx_irrep, irrep) in enumerate(irreps(psym))
            d = psym.character_table[idx_irrep, 1]
            @test irrep_dimension(psym, idx_irrep) == d
            for m in irrep
                @test size(m) == (d, d)
            end # for irrep
            for (idx_cc, cc) in enumerate(psym.conjugacy_classes)
                character = character_table(psym)[idx_irrep, idx_cc]
                for idx_elem in cc
                    @test isapprox(tr(irrep[idx_elem]), character; atol=Base.rtoldefault(Float64))
                end # for cc
            end # for enumerate(psym.conjugacy_classes)
        end # for enumerate(irreps(psym))
    end

    @testset "project" begin
        project(psym, [1 0 0; 0 1 0])
        @test_throws ArgumentError project(psym, [1 0; 0 1])
        @test_throws ArgumentError project(psym, [1 1 1; 1 1 1])
        @test_throws ArgumentError project(psym, [1 0 0;]) # projection not faithful
    end

    @testset "iscompatible" begin
        hc1 = HypercubicLattice([4 0; 0 4])
        hc2 = HypercubicLattice([4 0; 0 3])
        tsym1 = TranslationSymmetry(hc1)
        tsym2 = TranslationSymmetry(hc2)

        @test_throws DimensionMismatch iscompatible(hc1, psym)
        @test_throws DimensionMismatch iscompatible(tsym1, psym)

        psym_proj = project(psym, [1 0 0; 0 1 0])

        @test iscompatible(hc1, psym_proj)
        @test !iscompatible(hc2, psym_proj)
        @test iscompatible(tsym1, psym_proj)
        @test !iscompatible(tsym2, psym_proj)
    end


    @testset "two-band model" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
        addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

        psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

        @testset "lattice permutations" begin
            idx_C4 = 3
            @test element_name(psym, idx_C4) == "4<sup>+</sup><sub>001</sub>"
            @test findorbitalmap(unitcell, psym.matrix_representations[idx_C4]) == [(2, [0,0]), (1, [-1,0])]
            @test findorbitalmap(unitcell, psym)[idx_C4] == [(2, [0,0]), (1, [-1,0])]

            lattice = make_lattice(unitcell, [2 0; 0 2])
            tsym = TranslationSymmetry(lattice)

            perms = get_orbital_permutations(lattice, psym)
            @test length(perms) == length(psym.element_names)
            @test perms[1] == Permutation(1:8) # identity

            # |       |                 |       |
            # 6       8                 3       7
            # |       |         C4      |       |
            # . - 5 - . - 7 -   =>      . - 8 - . - 4 -
            # |       |                 |       |
            # 2       4                 1       5
            # |       |                 |       |
            # o - 1 - . - 3 -           o - 6 - . - 2 -
            @test perms[idx_C4] == Permutation([2,3,6,7,4,1,8,5])

        end # testset "lattice permutations"

        @testset "little group" begin
            lattice = make_lattice(unitcell, [2 0; 0 2])
            tsym = TranslationSymmetry(lattice)
            @test little_group_elements(tsym, 1, psym) == collect(1:8)
            @test little_group_elements(tsym, 2, psym) == [1,2,5,6]
            @test little_group_elements(tsym, 3, psym) == [1,2,5,6]
            @test little_group_elements(tsym, 4, psym) == collect(1:8)

            @test little_group(tsym, 1, psym) == little_group(tsym, psym, 1:8)
            @test little_group(tsym, 2, psym) == little_group(tsym, psym, [1,2,5,6])
            @test little_group(tsym, 3, psym) == little_group(tsym, psym, [1,2,5,6])
            @test little_group(tsym, 4, psym) == little_group(tsym, psym, 1:8)
        end

        @testset "little_symmetry" begin
            lattice = make_lattice(unitcell, [4 0; 0 4])
            tsym = TranslationSymmetry(lattice)
            for tsym_irrep in 1:num_irreps(tsym)
                psym_little = little_symmetry(tsym, tsym_irrep, psym)
                k = tsym.hypercube.coordinates[tsym_irrep]
                @test iscompatible(tsym, tsym_irrep, psym) == (k in [[0,0], [2,2]])
                @test iscompatible(tsym, tsym_irrep, psym_little)
                lg_matrep = psym.matrix_representations[little_group_elements(tsym, tsym_irrep, psym)]
                @test !isnothing(group_isomorphism(little_group(tsym, tsym_irrep, psym),
                                                                                      FiniteGroup(group_multiplication_table(lg_matrep))))
            end # for tsym_irrep
        end # testset little_symmetry
    end
end # @testset "PointSymmetry"
