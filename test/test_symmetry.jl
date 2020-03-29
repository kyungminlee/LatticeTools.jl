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
      @test tsym.hypercube.coordinates[idx_gen1] == [1, 0]
      @test tsym.hypercube.coordinates[idx_gen2] == [0, 1]
      # elements ordered according to the "generator" (i.e. orthogonal order)
      @test tsym.element_names == [
        "[0, 0]", "[1, 0]", "[2, 0]",
        "[0, 1]", "[1, 1]", "[2, 1]",
        "[0, 2]", "[1, 2]", "[2, 2]",
      ]
      @test tsym.conjugacy_classes == [(name=x, elements=[i]) for (i, x) in enumerate(tsym.element_names)]

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

  #=
  @testset "PointSymmetryBravaisRepresentation" begin
    psym = PointSymmetry(YAML.load_file("repr-19.yaml"))

    hypercube = HypercubicLattice([4 -2; 2 2])
    @test_throws ArgumentError PointSymmetryBravaisRepresentation(psym, hypercube, [1 0; 0 1])
    @test_throws ArgumentError PointSymmetryBravaisRepresentation(psym, hypercube, [1 1 1; 1 1 1])
    @test_throws ArgumentError PointSymmetryBravaisRepresentation(psym, HypercubicLattice([2 0; 0 1]), [1 0 0; 0 1 0])

    psbr = PointSymmetryBravaisRepresentation(psym, hypercube, [1 0 0; 0 1 0])

    n_elements = order(psbr.symmetry.group)
    for i in 1:n_elements, j in 1:n_elements
      k = findfirst(x->(x == psbr.permutations[i] * psbr.permutations[j]), psbr.permutations)
      @test psbr.symmetry.group.multiplication_table[i,j] == k
    end
  end
  =#


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
