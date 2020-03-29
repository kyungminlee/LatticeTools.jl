using Test
using LinearAlgebra
using YAML

@testset "Symmetry" begin
  @testset "TranslationSymmetry" begin

    @testset "orthogonal lattice" begin
      tsym = TranslationSymmetry([3 0; 0 3])

      @test length(tsym.generators) == 2
      idx_gen1 = tsym.generators[1]
      idx_gen2 = tsym.generators[2]

      @test group_order(tsym) == 9

      @testset "elements" begin
        @test tsym.hypercube.coordinates[idx_gen1] == [1, 0]
        @test tsym.hypercube.coordinates[idx_gen2] == [0, 1]
        # elements ordered according to the "generator" (i.e. orthogonal order)
        @test tsym.element_names == [
          "[0, 0]", "[1, 0]", "[2, 0]",
          "[0, 1]", "[1, 1]", "[2, 1]",
          "[0, 2]", "[1, 2]", "[2, 2]",
        ]
        @test length(element_names(tsym)) == 9
        @test element_names(tsym)[1] == element_name(tsym, 1)
        @test all(element_names(tsym)[i] == element_name(tsym, i) for i in 1:9)
      end

      @testset "coordinates" begin
        @test tsym.orthogonal_shape == [3, 3]
        @test tsym.orthogonal_coordinates == tsym.hypercube.coordinates
        for cc in tsym.hypercube.coordinates
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

      @testset "irreps" begin
        @test tsym.conjugacy_classes == [(name=x, elements=[i]) for (i, x) in enumerate(tsym.element_names)]
        @test size(character_table(tsym)) == (9, 9)
        @test length(irreps(tsym)) == 9
        @test all(irreps(tsym)[i] == irrep(tsym, i) for i in 1:9)
        @test num_irreps(tsym) == 9
        @test all(irrep_dimension(tsym, i) == 1 for i in 1:9)
      end
    end # @testset "orthogonal lattice" begin

    @testset "non-orthogonal lattice" begin
      tsym = TranslationSymmetry([4 0; 0 3])

      @test length(tsym.generators) == 1
      idx_gen = tsym.generators[1]
      @test tsym.hypercube.coordinates[idx_gen] == [1, 1]
      # elements ordered according to the "generator" (i.e. orthogonal order)
      @test tsym.element_names == [
        "[0, 0]", "[1, 1]", "[2, 2]", "[3, 0]",
        "[0, 1]", "[1, 2]", "[2, 0]", "[3, 1]",
        "[0, 2]", "[1, 0]", "[2, 1]", "[3, 2]",
      ]
      @test tsym.conjugacy_classes == [(name=x, elements=[i]) for (i, x) in enumerate(tsym.element_names)]

      @test tsym.orthogonal_shape == [12]
      @test tsym.orthogonal_coordinates == [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11]]
      for cc in tsym.hypercube.coordinates
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

  end # @testset "BravaisTranslationSymmetry" begin

  @testset "reduction" begin

    # Gamma point is always ok
    @test iscompatible([0,0], [3,3], [0,0])
    @test iscompatible([0,0], [3,3], [1,0])
    @test iscompatible([0,0], [3,3], [2,0])

    # non-zero momentum depends on what the identity translation is
    @test iscompatible([1,0], [3,3], [0,0]) # zero translation is always identity, so it's always fine
    @test !iscompatible([1,0], [3,3], [1,0]) # these two translations are not compatible
    @test !iscompatible([1,0], [3,3], [2,0]) #   with momentum [1,1]


    @test iscompatible([0,0], [3,3], [[0,0], [1,0], [2,0]])
    @test !iscompatible([1,0], [3,3], [[0,0], [1,0], [2,0]])
  end


  @testset "PointSymmetry" begin
    psym = PointSymmetryDatabase.get(19)
    @test psym.generators == [2,4]
    @test length(psym.conjugacy_classes) == 3 # three conjugacy classes
    @test size(psym.character_table) == (3,3)
    @test psym.character_table ≈ [1 1 1; 1 -1 1; 2 0 -1]
    @test all(let d = psym.character_table[idx_irrep, 1]
                size(m) == (d, d)
              end
                for (idx_irrep, irrep) in enumerate(irreps(psym))
                for m in irrep.matrices
              )
    ord_group = group_order(psym.group)
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

  end # @testset "PointSymmetry" begin

  @testset "PointSymmetry2" begin
    # D_4
    file_path = abspath(@__DIR__, "..", "data", "PointGroup3D", "PointGroup3D-12.yaml")
    data_yaml = YAML.load_file(file_path)

    let data_yaml_2 = deepcopy(data_yaml)
      data_yaml_2["Generators"] = [1]
      @test_throws ArgumentError PointSymmetry(data_yaml_2)
    end

    psym = PointSymmetry(data_yaml)
    @test group_order(psym) == 8
    let elnames = element_names(psym)
      for i in 1:8
        elnames[i] == element_name(psym, i)
      end
    end

    @test num_irreps(psym) == 5

    let χ = [1  1  1  1  1;
             1  1  1 -1 -1;
             1  1 -1 -1  1;
             1  1 -1  1 -1;
             2 -2  0  0  0]
      @test character_table(psym) ≈ χ
      @test size(character_table(psym)) == (5, 5)
    end

    let ir = irreps(psym)
      @test length(ir) == 5
      for i in 1:5
        @test ir[i] == irrep(psym, i)
      end
    end


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

  end # @testset "PointSymmetry"


  @testset "lattice permutation" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
    addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
    addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

    lattice = make_lattice(unitcell, [4 0; 0 4])

    tsym = TranslationSymmetry(lattice)

    perms = get_orbital_permutations(lattice)
    @test length(perms) == 16
    #@test length(perms[1].map) == 32
    @test perms[1].map == 1:32  # identity

    for (ir, r) in enumerate(lattice.hypercube.coordinates)
      @test perms[ir].map == [
        let
          ivec = [(i-1) % 4, (i-1) ÷ 4]
          jvec = [(x + 16) % 4 for x in (ivec .+ r)]
          j = jvec[1] + jvec[2]*4 + 1
          2 * (j-1) + iorb
        end
        for i in 1:16 for iorb in 1:2
      ]
    end

    @test_throws BoundsError get_irrep_iterator(lattice, tsym, 1, 2)
    let phases = [cis(-2π * i / 4) for j in 0:3 for i in 0:3]
      irrep_list1 = collect(get_irrep_iterator(lattice, tsym, 2, 1))
      irrep_list2 = collect(zip(perms, phases))
      @test length(irrep_list1) == length(irrep_list2)
      @test all(x[1] == y[1] for (x,y) in zip(irrep_list1, irrep_list2))
      @test all(x[2] ≈ y[2] for (x,y) in zip(irrep_list1, irrep_list2))
    end

    @test iscompatible(lattice, tsym, 1, [1,0]) # Γ point
    @test !iscompatible(lattice, tsym, 2, [1,0]) # Γ point
    @test !iscompatible(lattice, tsym, 2, [[0,0], [1,0]])


  end # testset lattice permutation



end # @testset "Symmetry" begin



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
#   @test Set(g.elements) == Set([t1^d1*t2^d2 for d1 in 0:2 for d2 in 0:1])
#   @test length(Set(g.elements)) == 2*3
#   @test g.fractional_momenta == [[0//3, 0//2], [1//3, 0//2], [2//3, 0//2],
#                                  [0//3, 1//2], [1//3, 1//2], [2//3, 1//2]]
#   χ = [cis(-2π * (k ⋅ t)) for k in g.fractional_momenta, t in g.translations]
#   @test isapprox(g.character_table, χ; atol=sqrt(eps(Float64)))
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
