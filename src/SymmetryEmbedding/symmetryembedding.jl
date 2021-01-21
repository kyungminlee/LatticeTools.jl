export AbstractSymmetryEmbedding

export SymmetryEmbedding

export embed
export element, elements
export symmetry
export generator_elements, generator_indices
export fractional_momentum

abstract type AbstractSymmetryEmbedding end

SymmetryOrEmbedding = Union{<:AbstractSymmetry, <:AbstractSymmetryEmbedding}

"""
    SymmetryEmbedding{SymmetryType<:AbstractSymmetry}

# Fields
- `lattice::Lattice`
- `symmetry::SymmetryType`
- `elements::Vector{SitePermutation}`
"""
struct SymmetryEmbedding{SymmetryType}<:AbstractSymmetryEmbedding
    lattice::Lattice
    symmetry::SymmetryType
    elements::Vector{SitePermutation}

    @doc """
        SymmetryEmbedding(lattice::Lattice, symmetry::SymmetryType)

    Construct a `SymmetryEmbedding`.
    """
    function SymmetryEmbedding(
        lattice::Lattice,
        symmetry::SymmetryType,
    ) where {SymmetryType<:AbstractSymmetry}
        if !iscompatible(lattice, symmetry)
            throw(ArgumentError("lattice and symmetry are not compatible"))
        end
        elems = [embed(lattice, elem) for elem in elements(symmetry)]
        if !allunique(elems)
            if isa(symmetry, PointSymmetry)
                reduced_elements = SitePermutation[]
                for elem in elems
                    if !in(elem, reduced_elements)
                        push!(reduced_elements, elem)
                    end
                end
                reduced_group = FiniteGroup(group_multiplication_table(reduced_elements))
                isomorphic_group_names = String[]
                for i in 1:32
                    target_symmetry = PointSymmetryDatabase.get(i)
                    if !isnothing(group_isomorphism(reduced_group, target_symmetry.group))
                        push!(isomorphic_group_names, target_symmetry.hermann_mauguin)
                    end
                end
                @warn "Lattice $(lattice.hypercube.shape_matrix) is too small for $(symmetry_name(symmetry)) (embedding not faithful).\n"*
                      "Reduced point group isomorphic to: $(join(isomorphic_group_names, ", "))" # COV_EXCL_LINE
            else
                @warn "Lattice $(lattice.hypercube.shape_matrix) is too small for $(symmetry_name(symmetry)) (embedding not faithful)." # COV_EXCL_LINE
            end
        end
        new{SymmetryType}(lattice, symmetry, elems)
    end
end

Base.eltype(::SymmetryEmbedding) = SitePermutation
Base.valtype(::SymmetryEmbedding) = SitePermutation
Base.eltype(::Type{<:SymmetryEmbedding}) = SitePermutation
Base.valtype(::Type{<:SymmetryEmbedding}) = SitePermutation

"""
    elements(symbed::SymmetryEmbedding)

Get elements of the symmetry embedding
"""
elements(symbed::SymmetryEmbedding) = symbed.elements

"""
    element(symbed::SymmetryEmbedding, g...)

Get the `g`th element of the symmetry embedding.
"""
element(symbed::SymmetryEmbedding, g...) = symbed.elements[g...]

"""
    symmetry(symbed::SymmetryEmbedding)

Return the underlying symmetry of the symmetry embedding
"""
symmetry(symbed::SymmetryEmbedding) = symbed.symmetry

"""
    generator_elements(symbed::SymmetryEmbedding)

Get the generators of the symmetry embedding, based on the generators of the underlying
symmetry. The embedding is allowed to be one-to-one, in which case the return value may
contain duplicates.
"""
generator_elements(symbed::SymmetryEmbedding) = element(symbed, generator_indices(symbed))


function Base.:(==)(lhs::SymmetryEmbedding{S}, rhs::SymmetryEmbedding{S}) where S
    return lhs.lattice == rhs.lattice && lhs.symmetry == rhs.symmetry && lhs.elements == rhs.elements
end


Base.iterate(sym::SymmetryEmbedding) = iterate(elements(sym))
Base.iterate(sym::SymmetryEmbedding, i) = iterate(elements(sym), i)


Base.length(sym::SymmetryEmbedding) = length(elements(sym))


for f in [
    :group_order,
    :group_multiplication_table,
    :generator_indices,
    :element_names,
    :element_name,
    :character_table,
    :irreps,
    :irrep,
    :irrep_dimension,
    :num_irreps, :numirreps, :irrepcount,
]
    eval(quote
        @doc """
            $($f)(::SymmetryEmbedding, args...)

        Calls `$($f)` for the underlying symmetry.
        See [`$($f)(::TranslationSymmetry)`](@ref), [`$($f)(::PointSymmetry)`](@ref).
        """
        ($f)(symbed::SymmetryEmbedding, args...) = ($f)(symmetry(symbed), args...)
    end)
end


"""
    fractional_momentum(symbed::SymmetryEmbedding{TranslationSymmetry}, args...)

Return fractional momentum(s) of the translation symmetry embedding.
"""
function fractional_momentum(symbed::SymmetryEmbedding{TranslationSymmetry}, args...)
    return fractional_momentum(symmetry(symbed), args...)
end


"""
    iscompatible(tsymbed, psymbed)

Check whether the point symmetry embedding `psymbed` is compatible with
the translation symmetry embedding `tsymbed`, i.e. whether they have
the same "lattice".
"""
function iscompatible(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    psymbed::SymmetryEmbedding{PointSymmetry},
)::Bool
    return tsymbed.lattice == psymbed.lattice
end


"""
    iscompatible(tsymbed, tsym_irrep_index, psymbed)

Check whether the point symmetry embedding `psymbed` is compatible with
the translation symmetry irrep defined by `tsym_irrep_index` and `symmetry(tsymbed)`.
In other words, the little group elements
"""
function iscompatible(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    tsym_irrep_index::Integer,
    psymbed::SymmetryEmbedding{PointSymmetry},
)::Bool
    !iscompatible(tsymbed, psymbed) && return false
    return little_group_elements(tsymbed, tsym_irrep_index, psymbed) == 1:group_order(psymbed)
end


"""
    little_group_elements(tsymbed, psymbed)
"""
function little_group_elements(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    psymbed::SymmetryEmbedding{PointSymmetry},
)
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    return little_group_elements(symmetry(tsymbed), symmetry(psymbed))
end


"""
    little_group_elements(tsymbed, tsym_irrep_index, psymbed)
"""
function little_group_elements(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    tsym_irrep_index::Integer,
    psymbed::SymmetryEmbedding{PointSymmetry},
)
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    return little_group_elements(symmetry(tsymbed), tsym_irrep_index, symmetry(psymbed))
end


"""
    little_symmetry(tsymbed, psymbed)
"""
function little_symmetry(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    psymbed::SymmetryEmbedding{PointSymmetry},
)
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    psym_little = little_symmetry(symmetry(tsymbed), symmetry(psymbed))
    return SymmetryEmbedding(psymbed.lattice, psym_little)
end


# export little_symmetry_strong

# WIP
# """
#     little_symmetry_strong(tsymbed, psymbed)
# """
# function little_symmetry_strong(tsymbed::SymmetryEmbedding{TranslationSymmetry},
#                                 psymbed::SymmetryEmbedding{PointSymmetry})
#     if !iscompatible(tsymbed, psymbed)
#         throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
#     end
#     psym_little = little_symmetry(symmetry(tsymbed), symmetry(psymbed))

#     little_element_indices = Int[]
#     little_element_embedding = Set{SitePermutation}()

#     for (i_elem, elem) in enumerate(elements(psym_little))
#         embed_elem = embed(tsymbed.lattice, elem)
#         if i_elem == 1
#             @assert iscompatible(symmetry(tsymbed), elem)
#             push!(little_element_indices, i_elem)
#             push!(little_element_embedding, embed_elem)
#         elseif (iscompatible(symmetry(tsymbed), elem) &&
#                 !isidentity(embed_elem) &&
#                 !in(embed_elem, little_element_embedding))
#             push!(little_element_indices, i_elem)
#             change = true
#             while change
#                 change = false
#                 push!(little_element_embedding, embed_elem)
#                 for elem1 in elements(tsymbed), elem2 in little_element_embedding
#                     elem3 = elem1 * elem2
#                     if !in(elem3, little_element_embedding)
#                         change = true
#                         push!(reduced_group_elements, elem3)
#                         break
#                     end
#                 end
#             end
#         end
#     end
#     @show little_element_indices
#     @show element_name(symmetry(psymbed), little_element_indices)
#     psym_little2 = little_symmetry(symmetry(tsymbed), symmetry(psymbed), little_element_indices)
#     return SymmetryEmbedding(psymbed.lattice, psym_little2)
# end




"""
    little_symmetry(tsymbed, tsym_irrep_index, psymbed)
"""
function little_symmetry(
    tsymbed::SymmetryEmbedding{TranslationSymmetry},
    tsym_irrep_index::Integer,
    psymbed::SymmetryEmbedding{PointSymmetry},
)
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    psym_little = little_symmetry(symmetry(tsymbed), tsym_irrep_index, symmetry(psymbed))
    return SymmetryEmbedding(psymbed.lattice, psym_little)
end


"""
    symmetry_name(arg::SymmetryEmbedding)

Return the name of the symmetry embedding.
"""
function symmetry_name(arg::SymmetryEmbedding)
    return "Embed[$(symmetry_name(symmetry(arg))) on $(arg.lattice.hypercube.shape_matrix)"*
        " with $(numsite(arg.lattice.unitcell)) sites]"
end


"""
    embed(lattice::Lattice, tsym::TranslationSymmetry)

Embed translation symmetry `tsym` to lattice.
"""
embed(lattice::Lattice, tsym::TranslationSymmetry) = SymmetryEmbedding(lattice, tsym)


"""
    embed(lattice::Lattice, psym::PointSymmetry)

Embed point symmetry `psym` to lattice.
"""
embed(lattice::Lattice, psym::PointSymmetry) = SymmetryEmbedding(lattice, psym)


"""
    translation_symmetry_embedding(lattice::Lattice)
"""
function translation_symmetry_embedding(lattice::Lattice)
    tsym = TranslationSymmetry(lattice)
    return embed(lattice, tsym)
end




# struct SymmorphicSpaceSymmetryEmbedding<:AbstractSymmetryEmbedding
#     lattice::Lattice
#     translation_symmetry::TranslationSymmetry
#     point_symmetry::PointSymmetry
#     elements::Matrix{SitePermutation}

#     function SymmorphicSpaceSymmetryEmbedding(
#                 lattice::Lattice,
#                 translation::TranslationSymmetry,
#                 point::PointSymmetry)
#         if !iscompatible(translation, point)
#             throw(ArgumentError("translation symmetry and point symmetry not compatible"))
#         end
#         tels = [embed(lattice, elem) for elem in translation.elements]
#         pels = [embed(lattice, elem) for elem in point.elements]

#         if !allunique(pels)
#             throw(ArgumentError("lattice too small for the point symmetry operation (not faithful)"))
#         end

#         elements = [pop * top for top in tels, pop in pels]
#         @assert size(elements) == (length(tels), length(pels))
#         new(lattice, translation, point, elements)
#     end
# end



# function embed(lattice::Lattice, tsym::TranslationSymmetry, psym::PointSymmetry)
#     SymmorphicSpaceSymmetryEmbedding(lattice, tsym, psym)
# end


# function symmetry_name(arg::SymmorphicSpaceSymmetryEmbedding)
#     name1 = symmetry_name(symmetry(arg.component1))
#     name2 = symmetry_name(symmetry(arg.component2))
#     return "Embed[$name1 ⋊ $name2 on $(arg.lattice.hypercube.scale_matrix)]"
# end
