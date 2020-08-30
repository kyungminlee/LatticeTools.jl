using Test
using LatticeTools

using LatticeTools: simplify_name

using LinearAlgebra
using YAML

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
        hermann_mauguin = "-1"
        schoenflies = "C<sub>i</sub>"
        psym = PointSymmetry(group, generators,
                             conjugacy_classes, character_table, irreps,
                             element_names, matrix_representations,
                             hermann_mauguin, schoenflies)

        let TBL = LatticeTools
            @test eltype(psym) == PointOperation{Int}
            @test valtype(psym) == PointOperation{Int}
            @test all(TBL.element(psym, i) == PointOperation(matrix_representations[i]) for i in 1:2)
            @test TBL.elements(psym) == PointOperation.(matrix_representations)
            @test all(TBL.element_name(psym, i) == element_names[i] for i in 1:2)
            @test TBL.element_names(psym) == element_names
            @test TBL.group(psym) == group
            @test TBL.group_order(psym) == group_order(group)
            @test TBL.group_multiplication_table(psym) == group.multiplication_table
            @test TBL.character_table(psym) == character_table
            @test TBL.irreps(psym) == irreps
            @test all(TBL.irrep(psym, i) == irreps[i] for i in 1:2)
            @test num_irreps(psym) == 2
            @test all(TBL.irrep_dimension(psym, i) == 1 for i in 1:2)
            @test generator_indices(psym) == [2]
            @test generator_elements(psym) == [PointOperation([-1 0; 0 -1])]
        end

        let generators = [1]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let conjugacy_classes = [[1], [1,2]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let character_table = [1 2 3; 4 5 6; 7 8 9]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [ones(ComplexF64, 1, 1), -ones(ComplexF64, 2, 2)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [2*ones(ComplexF64, 1, 1), -ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let irreps = [
                    [ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                    [ones(ComplexF64, 1, 1), 2*ones(ComplexF64, 1, 1)],
            ]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let element_names = ["1", "2", "3"]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let matrix_representations = [[1 0; 0 1], [-1 0 0; 0 -1 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let matrix_representations = [[1 0 0; 1 0 0], [-1 0 0; -1 0 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let matrix_representations = [[1 0 0; 1 0 0]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
        let matrix_representations = [[1 0 ; 0 1], [1 0; 0 1]]
            @test_throws ArgumentError PointSymmetry(group, generators,
                        conjugacy_classes, character_table, irreps,
                        element_names, matrix_representations, hermann_mauguin, schoenflies)
        end
    end

    # D_4 (422)
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

    @testset "symmetry_product" begin
        p = symmetry_product(psym)
        for x in elements(psym), y in elements(psym)
            z = x*y
            @test isa(z, PointOperation{Int})
            @test z.matrix == x.matrix * y.matrix
        end
    end

    @testset "iterate" begin
        @test [x for x in psym] == collect(elements(psym))
        @test all(x ∈ psym for x in psym)
        @test 100 ∉ psym
        @test PointOperation([1 2 0; -1 1 0; 0 0 1]) ∉ psym
        @test IdentityOperation(Int, 2) ∉ psym
        @test IdentityOperation(Int, 3) ∈ psym
        @test TranslationOperation([100, 100, 0]) ∉ psym
        @test TranslationOperation([0, 0]) ∉ psym
        @test TranslationOperation([0, 0, 0]) ∈ psym
        @test PointOperation([1 0 0; 0 1 0; 0 0 1]) ∈ psym
        @test PointOperation([0 1 0; 1 0 0; 0 0 -1]) ∈ psym
        @test SpaceOperation([0 1 0; 1 0 0; 0 0 -1], [0, 0, 0]) ∈ psym
        @test SpaceOperation([0 1 0; 1 0 0; 0 0 -1], [1, 0, 0]) ∉ psym
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
        hc1 = OrthoCube([4 0; 0 4])
        hc2 = OrthoCube([4 0; 0 3])
        tsym1 = TranslationSymmetry(hc1)
        tsym2 = TranslationSymmetry(hc2)

        @test_throws DimensionMismatch iscompatible(hc1, psym)
        @test_throws DimensionMismatch iscompatible(tsym1, psym)

        psym_proj = project(psym, [1 0 0; 0 1 0])

        @test iscompatible(hc1, psym_proj)
        @test !iscompatible(hc2, psym_proj)
        @test iscompatible(tsym1, psym_proj)
        @test !iscompatible(tsym2, psym_proj)
        @test little_symmetry(tsym1, psym_proj).hermann_mauguin == "422"
        @test little_symmetry(tsym2, psym_proj).hermann_mauguin == "222"
        @test little_symmetry(TranslationSymmetry([4 1; 0 3]), psym_proj).hermann_mauguin == "2"
        @test_throws ArgumentError little_symmetry(tsym2, 1, psym_proj) # when specifying irrep, tsym and psym have to be compatible
    end

    @testset "symmetry_name" begin
        let n = lowercase(symmetry_name(psym))
            @test occursin("point", n)
            @test occursin("422", n)
        end
    end


    @testset "two-band model" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
        addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
        addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

        psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

        @testset "lattice permutations" begin
            idx_C4 = 3
            @test element_name(psym, idx_C4) == "4<sup>+</sup><sub>001</sub>"
            @test findsitemap(unitcell, element(psym,idx_C4)) == [(2, [0,0]), (1, [-1,0])]
            @test findsitemap(unitcell, psym)[idx_C4] == [(2, [0,0]), (1, [-1,0])]

            lattice = make_lattice(unitcell, [2 0; 0 2])
            psymbed = embed(lattice, psym)
            tsym = TranslationSymmetry(lattice)
            tsymbed = embed(lattice, tsym)

            perms = [embed(lattice, op) for op in elements(psym)]
            @test length(perms) == length(psym.element_names)
            @test perms[1] == SitePermutation(1:8) # identity

            perms2 = SitePermutation[]
            for pop in elements(psym)
                push!(perms2, embed(lattice, pop))
                #push!(perms2, SitePermutation(get_site_permutation(lattice, mat, map)))
            end
            @test perms == perms2

            # |       |                 |       |
            # 6       8                 3       7
            # |       |         C4      |       |
            # . - 5 - . - 7 -   =>      . - 8 - . - 4 -
            # |       |                 |       |
            # 2       4                 1       5
            # |       |                 |       |
            # o - 1 - . - 3 -           o - 6 - . - 2 -
            @test perms[idx_C4] == SitePermutation([2,3,6,7,4,1,8,5])

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
            for LSYM in [little_symmetry, little_symmetry_iso]
                lattice = make_lattice(unitcell, [4 0; 0 4])
                tsym = TranslationSymmetry(lattice)
                for tsym_irrep in 1:num_irreps(tsym)
                    psym_little = LSYM(tsym, tsym_irrep, psym)
                    k = tsym.coordinates[tsym_irrep]
                    @test iscompatible(tsym, tsym_irrep, psym) == (k in [[0,0], [2,2]])
                    @test iscompatible(tsym, tsym_irrep, psym_little)
                    lg_matrep = psym.matrix_representations[little_group_elements(tsym, tsym_irrep, psym)]
                    @test !isnothing(group_isomorphism(little_group(tsym, tsym_irrep, psym),
                                                       FiniteGroup(group_multiplication_table(lg_matrep))
                                                       ))
                end # for tsym_irrep
            end
        end # testset little_symmetry
    end
end # @testset "PointSymmetry"


@testset "PointSymmetryDatabase" begin
    psym1 = PointSymmetryDatabase.get(13)
    psym2 = PointSymmetryDatabase.find("4mm")
    @test length(psym1) == length(psym2) == 8
    @test psym1.hermann_mauguin == "4mm"
    @test psym2.hermann_mauguin == "4mm"
    @test_throws ArgumentError PointSymmetryDatabase.get(999)
    @test isnothing(PointSymmetryDatabase.find("blah blah"))
    @test dimension(psym1) == dimension(psym2) == 3

    psym1p = PointSymmetryDatabase.get3d(13)
    psym2p = PointSymmetryDatabase.find3d("4mm")
    @test length(psym1p) == length(psym2p) == 8
    @test psym1p.hermann_mauguin == "4mm"
    @test psym2p.hermann_mauguin == "4mm"
    @test_throws ArgumentError PointSymmetryDatabase.get3d(999)
    @test isnothing(PointSymmetryDatabase.find3d("blah blah"))
    @test dimension(psym1p) == dimension(psym2p) == 3

    psym1p = PointSymmetryDatabase.get2d(6)
    psym2p = PointSymmetryDatabase.find2d("4mm")
    @test length(psym1p) == length(psym2p) == 8
    @test psym1p.hermann_mauguin == "4mm"
    @test psym2p.hermann_mauguin == "4mm"
    @test_throws ArgumentError PointSymmetryDatabase.get2d(999)
    @test isnothing(PointSymmetryDatabase.find2d("blah blah"))
    @test dimension(psym1p) == dimension(psym2p) == 2
end
