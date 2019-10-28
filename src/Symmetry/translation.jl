export TranslationGroup
export is_compatible
export minimal_translation_group
export get_generators
export translation_element

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
  return Permutation([hypercube.torus_wrap(r .+ displacement)[1] for r in hypercube.coordinates])
end

function groupmod(numer::Permutation, denominators ::AbstractVector{Permutation})
    min_perm = numer
    for denom in denominators
        g = numer * denom
        while g != numer
            if g < min_perm
                min_perm = g
            end
            g = g * denom
        end
    end
    return min_perm
end

function get_generators(hypercube ::HypercubicLattice)
  translations = [translation_element(hypercube, r) for r in hypercube.coordinates]
  indices = collect(1:length(hypercube.coordinates))
  generators = Permutation[]
  generator_indices = Int[]
  remaining_translations = Set(translations)
  function findmaxorder(indices ::Vector{Int}) ::Int # indices
    max_order = 1
    max_order_index = 1
    for i in indices
      g = translations[i]
      length(indices) % g.order != 0 && continue
      if g.order > max_order
        max_order = g.order
        max_order_index = i
      end
    end
    return max_order_index
  end

  while length(indices) > 1
    idx = findmaxorder(indices)
    @assert idx > 1
    push!(generators, translations[idx])
    push!(generator_indices, idx)
    remaining_translations = Set([groupmod(translations[i], generators) for i in indices])
    next_indices = Int[]
    for i in indices
      if translations[i] in remaining_translations && !(i in next_indices)
        push!(next_indices, i)
      end
    end
    indices = next_indices
  end
  return generator_indices
end


function minimal_translation_group(tentative_generators::AbstractArray{Permutation}; tol::Real=sqrt(eps(Float64)))
  all_elements = let
    shape = [g.order for g in tentative_generators]
    translations = vcat( collect( Iterators.product([0:g.order-1 for g in tentative_generators]...) )...)
    translations = [ [x...] for x in translations]
    Set([prod(gen^d for (gen, d) in zip(tentative_generators, dist)) for (ig, dist) in enumerate(translations)])
  end

  generators = Permutation[]
  for g in tentative_generators
    push!(generators, g)
    shape = [g.order for g in generators]
    translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
    translations = [ [x...] for x in translations]
    elements = Set([prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)])
    if all_elements == elements
      break
    end
  end
  return TranslationGroup(generators; tol=tol)
end

function minimal_translation_group(tentative_generators::Permutation...; tol::Real=sqrt(eps(Float64)))
  return minimal_translation_group([tentative_generators...])
end


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
