export AbstractSymmetryEmbedding

export SymmetryEmbedding
export SymmorphicSpaceSymmetryEmbedding

export embed
export element, elements
export symmetry
export generator_elements, generator_indices

abstract type AbstractSymmetryEmbedding<:AbstractSymmetry end

struct SymmetryEmbedding{SymmetryType}<:AbstractSymmetryEmbedding
    lattice::Lattice
    symmetry::SymmetryType
    elements::Vector{SitePermutation}

    function SymmetryEmbedding(lattice::Lattice, symmetry::SymmetryType) where {SymmetryType<:AbstractSymmetry}
        if !iscompatible(lattice, symmetry)
            throw(ArgumentError("lattice and symmetry are not compatible"))
        end
        elems = [embed(lattice, elem) for elem in elements(symmetry)]
        # @show symmetry
        # @show elems
        # @show hash.(elems)
        # @show allunique(elems)
        if !allunique(elems)
            throw(ArgumentError("lattice too small for the symmetry operation (not faithful)"))
        end
        new{SymmetryType}(lattice, symmetry, elems)
    end
end

elements(symbed::SymmetryEmbedding) = symbed.elements
element(symbed::SymmetryEmbedding, g) = symbed.elements[g]
symmetry(symbed::SymmetryEmbedding) = symbed.symmetry
generator_elements(symbed::SymmetryEmbedding) = element(symbed, generator_indices(symbed.symmetry))
generator_indices(symbed::SymmetryEmbedding) = generator_indices(symbed.symmetry)


for f in [:group_order,
          :group_multiplication_table,
          :element_names,
          :element_name,
          :character_table,
          :irreps,
          :irrep,
          :num_irreps,
          :irrep_dimension]
    eval(quote
        ($f)(symbed::SymmetryEmbedding, args...) = ($f)(symbed.symmetry, args...)
    end)
end


function iscompatible(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                      psymbed::SymmetryEmbedding{PointSymmetry})::Bool
    return tsymbed.lattice == psymbed.lattice
end


function iscompatible(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                      tsym_irrep_index::Integer,
                      psymbed::SymmetryEmbedding{PointSymmetry})::Bool
    # TODO: Check lattice?
    ! iscompatible(tsymbed, psymbed) && return false
    return little_group_elements(tsymbed, tsym_irrep_index, psymbed) == 1:group_order(psymbed)
end


function little_group_elements(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                               psymbed::SymmetryEmbedding{PointSymmetry})
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    return little_group_elements(symmetry(tsymbed), symmetry(psymbed))
end


function little_group_elements(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                               tsym_irrep_index::Integer,
                               psymbed::SymmetryEmbedding{PointSymmetry})
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    return little_group_elements(symmetry(tsymbed), tsym_irrep_index, symmetry(psymbed))
end


function little_symmetry(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                         psymbed::SymmetryEmbedding{PointSymmetry})
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    psym_little = little_symmetry(tsymbed.symmetry, psymbed.symmetry)
    return SymmetryEmbedding(psymbed.lattice, psym_little)
    # TODO: maybe the lattice is too small that psymbed elements become identity. Deal with those.
end


function little_symmetry(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                         tsym_irrep::Integer,
                         psymbed::SymmetryEmbedding{PointSymmetry})
    if !iscompatible(tsymbed, psymbed)
        throw(ArgumentError("translation and point symmetry-embeddings not compatible"))
    end
    psym_little = little_symmetry(symmetry(tsymbed), tsym_irrep, symmetry(psymbed))
    return SymmetryEmbedding(psymbed.lattice, psym_little)
end


function symmetry_name(arg::SymmetryEmbedding)
    return "Embed[$(symmetry_name(symmetry(arg))) on $(arg.lattice.hypercube.scale_matrix)]"
end


struct SymmorphicSpaceSymmetryEmbedding<:AbstractSymmetryEmbedding
    lattice::Lattice
    translation_symmetry::TranslationSymmetry
    point_symmetry::PointSymmetry
    elements::Matrix{SitePermutation}

    function SymmorphicSpaceSymmetryEmbedding(
                lattice::Lattice,
                translation::TranslationSymmetry,
                point::PointSymmetry)
        if !iscompatible(translation, point)
            throw(ArgumentError("translation symmetry and point symmetry not compatible"))
        end
        tels = [embed(lattice, elem) for elem in translation.elements]
        pels = [embed(lattice, elem) for elem in point.elements]

        if !allunique(pels)
            throw(ArgumentError("lattice too small for the point symmetry operation (not faithful)"))
        end

        elements = [pop * top for top in tels, pop in pels]
        @assert size(elements) == (length(tels), length(pels))
        new(lattice, translation, point, elements)
    end
end


embed(lattice::Lattice, tsym::TranslationSymmetry) = SymmetryEmbedding(lattice, tsym)
embed(lattice::Lattice, psym::PointSymmetry) = SymmetryEmbedding(lattice, psym)

function embed(lattice::Lattice, tsym::TranslationSymmetry, psym::PointSymmetry)
    SymmorphicSpaceSymmetryEmbedding(lattice, tsym, psym)
end


function symmetry_name(arg::SymmorphicSpaceSymmetryEmbedding)
    name1 = symmetry_name(symmetry(arg.component1))
    name2 = symmetry_name(symmetry(arg.component2))
    return "Embed[$name1 ⋊ $name2 on $(arg.lattice.hypercube.scale_matrix)]"
end
