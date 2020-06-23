using Test
using TightBindingLattice

using LinearAlgebra
using YAML


@testset "TranslationSymmetry" begin
    @testset "constructors" begin
        TranslationSymmetry(4*ones(Int, (1, 1)))
        TranslationSymmetry([4 0; 0 4])
        @test_throws ErrorException TranslationSymmetry([4 0 0; 0 4 0; 0 0 4])  # <- Temporary
        @test_throws ArgumentError TranslationSymmetry(OrthoCube([2 0; 0 2]), [2 0; 0 1])
    end

    @testset "properties" begin
        tsym = TranslationSymmetry([3 0; 0 1])
        mtab = [1 2 3; 2 3 1; 3 1 2]
        @test eltype(tsym) == TranslationOperation{Int}
        @test valtype(tsym) == TranslationOperation{Int}
        @test group(tsym) == FiniteGroup(mtab)
        @test group_order(tsym) == 3
        @test group_order(tsym, 1) == 1
        @test group_order(tsym, 2) == 3
        @test group_order(tsym, 3) == 3
        @test group_multiplication_table(tsym) == mtab
        @test element(tsym, 1) == TranslationOperation([0,0])
        @test element(tsym, 2) == TranslationOperation([1,0])
        @test element(tsym, 3) == TranslationOperation([2,0])
        @test elements(tsym) == TranslationOperation.([[0,0], [1,0], [2,0]])
        @test [element_name(tsym, i) for i in 1:3] == element_names(tsym)
        let ω = cis(2π/3)
            @test character_table(tsym) ≈ [1 1 1; 1 ω^2 ω; 1 ω ω^2]
        end
        @test length(irreps(tsym)) == 3
        @test num_irreps(tsym) == 3
        @test irrep(tsym, 1) == irreps(tsym)[1]
        @test irrep(tsym, 2) == irreps(tsym)[2]
        @test irrep(tsym, 3) == irreps(tsym)[3]
        @test all(irrep_dimension(tsym, i) == 1 for i in 1:3)
        @test generator_indices(tsym) == [2, 1]
        @test generator_elements(tsym) == [TranslationOperation([1,0]), TranslationOperation([0,0])]

        let n = lowercase(symmetry_name(tsym))
            @test occursin("translation", n)
            @test occursin("3", n)
            @test occursin("1", n)
            @test occursin("0", n)
        end
    end

    @testset "orthogonal lattice" begin
        tsym = TranslationSymmetry([3 0; 0 3])
        @test isabelian(tsym.group)

        @test length(tsym.generators) == 2
        idx_gen1 = tsym.generators[1]
        idx_gen2 = tsym.generators[2]

        @test group_order(tsym) == 9

        @testset "generators" begin
            # Generators completely generate the group
            @test generate_subgroup(tsym.group, tsym.generators) == BitSet(1:group_order(tsym))
            @test prod(tsym.group.period_lengths[x] for x in tsym.generators) == group_order(tsym) # since abelian group
        end

        @testset "elements" begin
            @test tsym.coordinates[idx_gen1] == [1, 0]
            @test tsym.coordinates[idx_gen2] == [0, 1]
            # elements ordered according to the "generator" (i.e. orthogonal order)
            @test tsym.element_names == [
                "[0, 0]", "[1, 0]", "[2, 0]",
                "[0, 1]", "[1, 1]", "[2, 1]",
                "[0, 2]", "[1, 2]", "[2, 2]",
            ]
            @test length(element_names(tsym)) == 9
            for i in 1:9
                @test element_names(tsym)[i] == element_name(tsym, i)
            end
        end

        @testset "multiplication table" begin
            mtab = zeros(Int, (9, 9))
            lookup = Dict(r => i for (i, r) in enumerate(tsym.coordinates))
            for (i, ri) in enumerate(tsym.coordinates)
                for (j, rj) in enumerate(tsym.coordinates)
                    _, ek = tsym.orthocube.wrap(ri + rj)
                    k = lookup[ek]
                    mtab[i,j] = k
                end # for j
            end # for i
            @test mtab == group_multiplication_table(tsym)
            @test mtab == group_multiplication_table(tsym.group)
        end

        @testset "coordinates" begin
            @test tsym.orthogonal_shape == [3, 3]
            @test tsym.orthogonal_coordinates == tsym.coordinates
            for cc in tsym.coordinates
                oc = tsym.coordinate_to_orthogonal_map[cc]
                cc2 = tsym.orthogonal_to_coordinate_map[oc]
                @test cc == cc2
            end

            for oc in tsym.orthogonal_coordinates
                cc = tsym.orthogonal_to_coordinate_map[oc]
                oc2 = tsym.coordinate_to_orthogonal_map[cc]
                @test oc == oc2
            end
        end

        @testset "symmetry_product" begin
            pro = symmetry_product(tsym)
            @test all( pro(x,y).displacement == tsym.orthocube.wrap(x.displacement + y.displacement)[2]
                       for x in tsym, y in tsym)
        end

        @testset "iterate/membership" begin
            @test [x for x in tsym] == collect(elements(tsym))
            @test all(x ∈ tsym for x in tsym)
            @test 100 ∉ tsym
            @test IdentityOperation(Int, 2) ∈ tsym
            @test IdentityOperation(Int, 3) ∉ tsym
            @test TranslationOperation([100, 100]) ∈ tsym
            @test TranslationOperation([100, 100, 100]) ∉ tsym
            @test PointOperation([1 0; 0 1]) ∈ tsym
            @test PointOperation([0 1; 1 0]) ∉ tsym
            @test SpaceOperation([1 0; 0 1], [1, 0]) ∈ tsym
            @test SpaceOperation([0 1; 1 0], [1, 0]) ∉ tsym
        end

        @testset "irreps" begin
            @test tsym.conjugacy_classes == [[i] for (i, x) in enumerate(tsym.element_names)]
            @test size(character_table(tsym)) == (9, 9)
            @test length(irreps(tsym)) == 9
            @test num_irreps(tsym) == 9
            for i in 1:9
                @test irreps(tsym)[i] == irrep(tsym, i)
                @test irrep_dimension(tsym, i) == 1
            end
        end
    end # @testset "orthogonal lattice" begin

    @testset "non-orthogonal lattice" begin
        tsym = TranslationSymmetry([4 0; 0 3])

        @test length(tsym.generators) == 2
        # idx_gen = tsym.generators[1]
        @test tsym.coordinates[tsym.generators[1]] == [1, 0]
        @test tsym.coordinates[tsym.generators[2]] == [0, 1]
        # elements ordered according to the "generator" (i.e. orthogonal order)
        @test element_names(tsym) == [
            "[0, 0]", "[1, 0]", "[2, 0]", "[3, 0]",
            "[0, 1]", "[1, 1]", "[2, 1]", "[3, 1]",
            "[0, 2]", "[1, 2]", "[2, 2]", "[3, 2]",
        ]
        @test tsym.conjugacy_classes == [[i] for (i, x) in enumerate(tsym.element_names)]

        @test tsym.orthogonal_shape == [4,3]
        @test tsym.orthogonal_coordinates == [[0,0], [1,0], [2,0], [3,0], [0,1], [1,1], [2,1], [3,1], [0,2], [1,2], [2,2], [3,2]]
        for cc in tsym.coordinates
            oc = tsym.coordinate_to_orthogonal_map[cc]
            cc2 = tsym.orthogonal_to_coordinate_map[oc]
            @test cc == cc2
        end

        for oc in tsym.orthogonal_coordinates
            cc = tsym.orthogonal_to_coordinate_map[oc]
            oc2 = tsym.coordinate_to_orthogonal_map[cc]
            @test oc == oc2
        end
    end  # @testset "non-orthogonal lattice" begin

    @testset "lattice permutation" begin
        unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
        addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
        lattice = make_lattice(unitcell, [4 0; 0 4])
        tsym = TranslationSymmetry(lattice)

        let
            lattice_failure = make_lattice(unitcell, [4 0; 0 3])
            #@test_throws ArgumentError get_site_permutations(lattice_failure, tsym)
            @test_throws ArgumentError embed(lattice_failure, tsym)
        end
        perms = [embed(lattice, op) for op in tsym.elements]
        #tsymbed = embed(lattice, tsym)
        @test length(perms) == 16
        @test perms[1].permutation.map == 1:32  # identity

        for (ir, r) in enumerate(lattice.bravais_coordinates)
            @test perms[ir].permutation.map == [
                let
                    ivec = [(i-1) % 4, (i-1) ÷ 4]
                    jvec = [(x + 16) % 4 for x in (ivec .+ r)]
                    j = jvec[1] + jvec[2]*4 + 1
                    2 * (j-1) + iorb
                end
                for i in 1:16 for iorb in 1:2
            ]
        end

        tsymbed = embed(lattice, tsym)
        @test_throws ArgumentError get_irrep_iterator(IrrepComponent(tsymbed, 1, 2))
        let phases = [cis(-2π * i / 4) for j in 0:3 for i in 0:3]
            tsic = IrrepComponent(tsymbed, 2, 1)
            irrep_list1 = collect(get_irrep_iterator(tsic))
            irrep_list2 = collect(zip(perms, phases))
            @test length(irrep_list1) == length(irrep_list2)
            for (x,y) in zip(irrep_list1, irrep_list2)
                @test x[1] == y[1]
                @test x[2] ≈ y[2]
            end
        end

        @test  isbragg(tsym.fractional_momenta[1], [1,0]) # Γ point
        @test !isbragg(tsym.fractional_momenta[2], [1,0]) # Γ point
        @test !isbragg(tsym.fractional_momenta[2], [[0,0], [1,0]])
    end # testset lattice permutation

    @testset "reduction" begin
        @test_throws DimensionMismatch isbragg([4, 0], [2,0,0], [2,0])
        @test_throws DimensionMismatch isbragg([4, 4], [2,0], [2,0,0])

        @test_throws ArgumentError isbragg([4, 0], [2,0], [2,0])
        @test_throws ArgumentError isbragg([-2, 2], [2,0], [2,0])

        tsym = TranslationSymmetry([3 0; 0 3])
        # Gamma point is always ok
        @test  isbragg([3,3], [0,0], [0,0])
        @test  isbragg([3,3], [0,0], [1,0])
        @test  isbragg([3,3], [0,0], [2,0])

        @test  isbragg(tsym, 1, TranslationOperation([0,0]))
        @test  isbragg(tsym, 1, TranslationOperation([1,0]))
        @test  isbragg(tsym, 1, TranslationOperation([2,0]))

        # non-zero momentum depends on what the identity translation is
        @test  isbragg([3,3], [1,0], [0,0]) # zero translation is always identity, so it's always fine
        @test !isbragg([3,3], [1,0], [1,0]) # these two translations are not compatible
        @test !isbragg([3,3], [1,0], [2,0]) #   with momentum [1,1]

        @test  isbragg(tsym, 2, TranslationOperation([0,0]))
        @test !isbragg(tsym, 2, TranslationOperation([1,0]))
        @test !isbragg(tsym, 2, TranslationOperation([2,0]))

        @test  isbragg([3,3], [0,0], [[0,0], [1,0], [2,0]])
        @test !isbragg([3,3], [1,0], [[0,0], [1,0], [2,0]])

        @test  isbragg(tsym, 1, TranslationOperation.([[0,0], [1,0], [2,0]]))
        @test !isbragg(tsym, 2, TranslationOperation.([[0,0], [1,0], [2,0]]))

        @test_throws DimensionMismatch isbragg([1//2, 2//4, 0//2], [2,0])

        # Gamma point is always ok
        @test  isbragg([0,0] .// [3,3], [0,0])
        @test  isbragg([0,0] .// [3,3], [1,0])
        @test  isbragg([0,0] .// [3,3], [2,0])

        # non-zero momentum depends on what the identity translation is
        @test  isbragg([1,0] .// [3,3], [0,0]) # zero translation is always identity, so it's always fine
        @test !isbragg([1,0] .// [3,3], [1,0]) # these two translations are not compatible
        @test !isbragg([1,0] .// [3,3], [2,0]) #   with momentum [1,1]

        @test  isbragg([0,0] .// [3,3], [[0,0], [1,0], [2,0]])
        @test !isbragg([1,0] .// [3,3], [[0,0], [1,0], [2,0]])


        tsym = TranslationSymmetry([4 0; 0 4])
        @test  isbragg(tsym, 11, TranslationOperation([2,0]))
        @test !isbragg(tsym, 10, TranslationOperation([2,0]))
        @test  isbragg(tsym, 11, TranslationOperation.([[0,0], [2,0], [0,2], [2,2]]))

        tsym = TranslationSymmetry([6 2; -2 6])

        # TODO more tests for isbragg

        for irrep_index in 1:num_irreps(tsym)
            k = tsym.fractional_momenta[irrep_index]
            k_ortho = tsym.orthogonal_coordinates[irrep_index] .// tsym.orthogonal_shape
            for (t, t_ortho) in zip(elements(tsym), tsym.orthogonal_coordinates)
                @test mod(k⋅t.displacement, 1) == mod(k_ortho⋅t_ortho, 1)
            end
        end
    end
end # @testset "TranslationSymmetry" begin



#
# @testset "translation" begin
#   @testset "constructor exceptions" begin
#     @test_throws ArgumentError TranslationGroup([Permutation([2,3,4,1]), Permutation([3,4,1,2])])
#     @test_throws ArgumentError TranslationGroup([Permutation([2,1,3,4]), Permutation([1,3,2,4])])
#   end
#
#   t1 = Permutation([2,3,1, 5,6,4])
#   t2 = Permutation([4,5,6, 1,2,3])
#   g = TranslationGroup([t1, t2])
#
#   @test g.generators == [t1, t2]
#   @test g.translations == [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1]]
#   @test Set(g) == Set([t1^d1*t2^d2 for d1 in 0:2 for d2 in 0:1])
#   @test length(Set(g)) == 2*3
#   @test g.fractional_momenta == [[0//3, 0//2], [1//3, 0//2], [2//3, 0//2],
#                                  [0//3, 1//2], [1//3, 1//2], [2//3, 1//2]]
#   χ = [cis(-2π * (k ⋅ t)) for k in g.fractional_momenta, t in g.translations]
#   @test isapprox(g.character_table, χ; atol=Base.rtoldefault(Float64))
#
#   @test is_compatible([0//1, 0//1], [0,0])
#   @test !is_compatible([0//1, 1//2], [0,1])
#
#   @test is_compatible([0//1, 0//1], [[0,0]])
#   @test !is_compatible([0//1, 1//2], [[0,0], [0,1]])
# end
#
#
# @testset "non-orthogonal" begin
#   #=
#   . . . o .
#   o . . . .
#   . . . . o
#   . o . . .
#
#   2 3 4 1 5
#   1|5 6 7|2
#   7|2 3 4|1
#   4|1|5 6 7
#
#
#   =#
#   t1 = Permutation([5,3,4,1,6,7,2])
#   t2 = Permutation([2,5,6,7,3,4,1])
#   generators = [t1, t2]
#   #translation_group = TranslationGroup(t1, t2)
#   full_shape = [g.order for g in generators]
#   full_translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
#   full_translations = [ [x...] for x in full_translations]
#
#   elemtrans = Dict{Permutation, Vector{Int}}()
#   for (ig, dist) in enumerate(full_translations)
#     g = prod(gen^d for (gen, d) in zip(generators, dist))
#     if !haskey(elemtrans, g)
#       elemtrans[g] = [dist...]
#     end
#   end
#   elements = Permutation[]
#   translations = Vector{Int}[]
#   for (e, t) in elemtrans
#     push!(elements, e)
#     push!(translations, t)
#   end
#   idx = sortperm(translations)
#   translations = translations[idx]
#   elements = elements[idx]
# end
