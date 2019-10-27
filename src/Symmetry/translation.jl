export TranslationGroup
export is_compatible
export minimal_translation_group

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

    shape = [g.cycle_length for g in generators]
    translations = vcat( collect( Iterators.product([0:g.cycle_length-1 for g in generators]...) )...)
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

function minimal_translation_group(tentative_generators::AbstractArray{Permutation}; tol::Real=sqrt(eps(Float64)))
  all_elements = let
    shape = [g.cycle_length for g in tentative_generators]
    translations = vcat( collect( Iterators.product([0:g.cycle_length-1 for g in tentative_generators]...) )...)
    translations = [ [x...] for x in translations]
    Set([prod(gen^d for (gen, d) in zip(tentative_generators, dist)) for (ig, dist) in enumerate(translations)])
  end

  generators = Permutation[]
  for g in tentative_generators
    push!(generators, g)
    shape = [g.cycle_length for g in generators]
    translations = vcat( collect( Iterators.product([0:g.cycle_length-1 for g in generators]...) )...)
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
