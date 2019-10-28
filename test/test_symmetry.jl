using Test

using LinearAlgebra

#include("pauli_matrix.jl")

@testset "Permutation" begin
  @test_throws ArgumentError Permutation([1,2,4])
  @test_throws OverflowError Permutation([mod(x, 4096)+1 for x in 1:4096])
  p0 = Permutation([1,2,3,4])
  p1 = Permutation([2,3,4,1])
  p2 = Permutation([3,4,1,2])
  p3 = Permutation([4,1,2,3])

  @test p1 * p2 == p3
  @test p1 != p3
  @test p1^0 == p0
  @test p1^1 == p1
  @test p1^2 == p2
  @test p1^3 == p3

  @test p0.order == 1
  @test p1.order == 4
  @test p2.order == 2
  @test p3.order == 4

  @test_throws ArgumentError Permutation([1,2,3,4]) * Permutation([1,2,3,4,5])
  @test hash(Permutation(Int[1,2,3,4])) == hash(Int[1,2,3,4])
end

@testset "translation" begin
  @testset "constructor exceptions" begin
    @test_throws ArgumentError TranslationGroup([Permutation([2,3,4,1]), Permutation([3,4,1,2])])
    @test_throws ArgumentError TranslationGroup([Permutation([2,1,3,4]), Permutation([1,3,2,4])])
  end

  t1 = Permutation([2,3,1, 5,6,4])
  t2 = Permutation([4,5,6, 1,2,3])
  g = TranslationGroup([t1, t2])

  @test g.generators == [t1, t2]
  @test g.translations == [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1]]
  @test Set(g.elements) == Set([t1^d1*t2^d2 for d1 in 0:2 for d2 in 0:1])
  @test length(Set(g.elements)) == 2*3
  @test g.fractional_momenta == [[0//3, 0//2], [1//3, 0//2], [2//3, 0//2],
                                 [0//3, 1//2], [1//3, 1//2], [2//3, 1//2]]
  χ = [cis(2π * (k ⋅ t)) for k in g.fractional_momenta, t in g.translations]
  @test isapprox(g.character_table, χ; atol=sqrt(eps(Float64)))

  @test is_compatible([0//1, 0//1], [0,0])
  @test !is_compatible([0//1, 1//2], [0,1])

  @test is_compatible([0//1, 0//1], [[0,0]])
  @test !is_compatible([0//1, 1//2], [[0,0], [0,1]])
end


@testset "non-orthogonal" begin
  #=
  . . . o .
  o . . . .
  . . . . o
  . o . . .

  2 3 4 1 5
  1|5 6 7|2
  7|2 3 4|1
  4|1|5 6 7


  =#
  t1 = Permutation([5,3,4,1,6,7,2])
  t2 = Permutation([2,5,6,7,3,4,1])
  generators = [t1, t2]
  #translation_group = TranslationGroup(t1, t2)
  full_shape = [g.order for g in generators]
  full_translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
  full_translations = [ [x...] for x in full_translations]

  elemtrans = Dict{Permutation, Vector{Int}}()
  for (ig, dist) in enumerate(full_translations)
    g = prod(gen^d for (gen, d) in zip(generators, dist))
    if !haskey(elemtrans, g)
      elemtrans[g] = [dist...]
    end
  end
  elements = Permutation[]
  translations = Vector{Int}[]
  for (e, t) in elemtrans
    push!(elements, e)
    push!(translations, t)
  end
  idx = sortperm(translations)
  translations = translations[idx]
  elements = elements[idx]
end
