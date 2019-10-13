export TranslationGroup
export is_compatible

struct TranslationGroup <: AbstractSymmetryGroup
  generators ::Vector{Permutation}
  
  translations ::Vector{Vector{Int}}
  elements ::Vector{Permutation}
  fractional_momenta ::Vector{Vector{Rational}}
  
  conjugacy_classes ::Vector{Int}  # conjugacy class of element
  character_table ::Array{ComplexF64, 2}

  function TranslationGroup(generators::AbstractArray{Permutation})
    @assert all(g1 * g2 == g2 * g1 for g1 in generators, g2 in generators) "non-commuting set of generators"

    shape = [g.cycle_length for g in generators]
    translations = vcat( collect( Iterators.product([0:g.cycle_length-1 for g in generators]...) )...)
    translations = [ [x...] for x in translations]
    elements = [prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)]

    @assert length(Set(elements)) == length(elements) "elements not unique (generators not orthogonal)"

    momentum(sub) = [x//d for (x, d) in zip(sub, shape)]
    fractional_momenta = [momentum(sub) for sub in translations]

    conjugacy_classes = collect(1:length(elements))
    character_table = ComplexF64[cis(dot(float.(kf) .* 2π, t))
                                 for kf in fractional_momenta, t in translations]

    return new(generators, translations, elements, fractional_momenta, conjugacy_classes, character_table)
  end
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
    fractional_momentum ::AbstractVector{Rational},
    identity_translation ::AbstractVector{<:Integer}
    )
  value = sum( i * j for (i,j) in zip(fractional_momentum, identity_translation))
  return mod(value, 1) == 0
end


function is_compatible(
    fractional_momentum ::AbstractVector{Rational},
    identity_translations ::AbstractVector{<:AbstractVector{<:Integer}}
  )
  return all(is_compatible(fractional_momentum, t) for t in identity_translations)
end

