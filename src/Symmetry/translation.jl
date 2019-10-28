export TranslationGroup
export is_compatible
#export minimal_translation_group
export get_generators
export translation_element
export translation_symmetry_group

struct TranslationGroup <: AbstractSymmetryGroup
  generators ::Vector{Permutation}

  translations ::Vector{Vector{Int}}
  elements ::Vector{Permutation}
  fractional_momenta ::Vector{Vector{Rational}}

  conjugacy_classes ::Vector{Int}  # conjugacy class of element
  character_table ::Array{ComplexF64, 2}

  element_irreps ::Array{Array{ComplexF64, 2}, 2}

  function TranslationGroup(p::Permutation...; tol::Real=sqrt(eps(Float64)))
    return TranslationGroup(Permutation[p...])
  end

  function TranslationGroup(generators::AbstractArray{Permutation}; tol::Real=sqrt(eps(Float64)))
    if ! all(g1 * g2 == g2 * g1 for g1 in generators, g2 in generators)
      throw(ArgumentError("non-commuting set of generators"))
    end

    shape = [g.order for g in generators]
    translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
    translations = [ [x...] for x in translations]
    elements = [prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)]

    if length(Set(elements)) != length(elements)
      throw(ArgumentError("elements not unique (generators not orthogonal)"))
    end

    momentum(sub) = [x//d for (x, d) in zip(sub, shape)]
    fractional_momenta = [momentum(sub) for sub in translations]

    conjugacy_classes = collect(1:length(elements))
    character_table = ComplexF64[cis(dot(float.(kf) .* 2π, t))
                                 for kf in fractional_momenta, t in translations]

    character_table_r = real.(character_table)
    character_table_i = imag.(character_table)

    character_table_r[abs.(character_table_r) .< tol] .= 0.0
    character_table_i[abs.(character_table_i) .< tol] .= 0.0

    character_table = character_table_r + im * character_table_i

    element_irreps = [ c * ones(ComplexF64, 1, 1) for c in character_table ]

    return new(generators, translations, elements, fractional_momenta, conjugacy_classes, character_table, element_irreps)
  end
end


function translation_element(hypercube ::HypercubicLattice,
                             displacement ::AbstractVector{<:Integer})
  return Permutation([hypercube.torus_wrap(r .+ displacement)[2] for r in hypercube.coordinates])
end


function get_generators(hypercube ::HypercubicLattice)
  translations = [translation_element(hypercube, r) for r in hypercube.coordinates]
  group_order = length(translations)
  sorted_indices = sortperm(1:group_order; by=(i->translations[i].order), rev=true, alg=Base.Sort.DEFAULT_STABLE)
  function factorize(generator_indices ::AbstractVector{<:Integer}, # Return value if successful
                     translation_groups ::AbstractVector{<:AbstractSet{Permutation}},
                     size_product ::Integer,
                     sorted_indices_begin ::Integer)    # ORDER BY "order" DESC
    # checking identity is enough
    length(translation_groups[1]) == group_order && return generator_indices
    for i in sorted_indices_begin:group_order
      idx = sorted_indices[i]
      g = translations[idx]
      if group_order % (size_product * g.order) == 0
        new_translation_groups = [tg * g for tg in translation_groups]
        push!(generator_indices, idx)
        ret = factorize(generator_indices,
                        new_translation_groups,
                        size_product * g.order,
                        i+1)
        ret !== nothing && return ret
        pop!(generator_indices)
      end
    end
    return nothing
  end
  generator_indices = Int[]
  sizehint!(generator_indices, group_order) # upper limit
  initial_translation_groups = [generate_group(g) for g in translations]
  return factorize(generator_indices, initial_translation_groups, 1, 1)
end


export translation_symmetry_group
function translation_symmetry_group(hypercube ::HypercubicLattice)
  generator_indices = get_generators(hypercube)
  generators = [translation_element(hypercube, hypercube.coordinates[i]) for i in generator_indices]
  return TranslationGroup(generators)
end



#
# function minimal_translation_group(tentative_generators::AbstractArray{Permutation}; tol::Real=sqrt(eps(Float64)))
#   all_elements = let
#     shape = [g.order for g in tentative_generators]
#     translations = vcat( collect( Iterators.product([0:g.order-1 for g in tentative_generators]...) )...)
#     translations = [ [x...] for x in translations]
#     Set([prod(gen^d for (gen, d) in zip(tentative_generators, dist)) for (ig, dist) in enumerate(translations)])
#   end
#
#   generators = Permutation[]
#   for g in tentative_generators
#     push!(generators, g)
#     shape = [g.order for g in generators]
#     translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
#     translations = [ [x...] for x in translations]
#     elements = Set([prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)])
#     if all_elements == elements
#       break
#     end
#   end
#   return TranslationGroup(generators; tol=tol)
# end
#
# function minimal_translation_group(tentative_generators::Permutation...; tol::Real=sqrt(eps(Float64)))
#   return minimal_translation_group([tentative_generators...])
# end
#

"""
    is_compatible

Check whether the fractional momentum ([0, 1)ᴺ) compatible with the identity translation.
i.e. k¹ R¹ + k² R² + ... + kᴺ Rᴺ = 0 (mod 1)

# Arguments
- `fractional_momentum ::AbstractVector{Rational}` : k
- `identity_translation ::AbstractVector{<:Integer}` : R
"""
function is_compatible(
    fractional_momentum ::AbstractVector{<:Rational},
    identity_translation ::AbstractVector{<:Integer}
    )
  value = sum( i * j for (i,j) in zip(fractional_momentum, identity_translation))
  return mod(value, 1) == 0
end


function is_compatible(
    fractional_momentum ::AbstractVector{<:Rational},
    identity_translations ::AbstractVector{<:AbstractVector{<:Integer}}
  )
  return all(is_compatible(fractional_momentum, t) for t in identity_translations)
end
