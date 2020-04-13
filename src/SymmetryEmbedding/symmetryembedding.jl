export AbstractSymmetryEmbedding

export SymmetryEmbedding
export SymmorphicSpaceSymmetryEmbedding

export embed
export element, elements
export symmetry

abstract type AbstractSymmetryEmbedding <:AbstractSymmetry end

struct SymmetryEmbedding{SymmetryType} <:AbstractSymmetryEmbedding
    lattice::Lattice
    symmetry::SymmetryType
    elements::Vector{SitePermutation}

    function SymmetryEmbedding(lattice::Lattice, symmetry::SymmetryType) where {SymmetryType<:AbstractSymmetry}
        if !iscompatible(lattice, symmetry)
            throw(ArgumentError("lattice and symmetry are not compatible"))
        end
        elements = [embed(lattice, elem) for elem in symmetry.elements]
        new{SymmetryType}(lattice, symmetry, elements)
    end
end

elements(symbed::SymmetryEmbedding) = symbed.elements
element(symbed::SymmetryEmbedding, g) = symbed.elements[g]
symmetry(symbed::SymmetryEmbedding) = symbed.symmetry

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

# struct TranslationSymmetryEmbedding<:AbstractSymmetryEmbedding
#     lattice::Lattice
#     symmetry::TranslationSymmetry
#     elements::Vector{SitePermutation}

#     function TranslationSymmetryEmbedding(lattice::Lattice, symmetry::TranslationSymmetry)
#         if (lattice.hypercube != symmetry.hypercube)
#             throw(ArgumentError("lattice and translation symmetry not compatible"))
#         end
#         elements = [embed(lattice, elem) for elem in symmetry.elements]
#         new(lattice, symmetry, elements)
#     end
# end


# elements(symbed::TranslationSymmetryEmbedding) = symbed.elements
# element(symbed::TranslationSymmetryEmbedding, g) = symbed.elements[g]


# struct PointSymmetryEmbedding<:AbstractSymmetryEmbedding
#     lattice::Lattice
#     symmetry::PointSymmetry
#     elements::Vector{SitePermutation}

#     function PointSymmetryEmbedding(lattice::Lattice, symmetry::PointSymmetry)
#         if !iscompatible(lattice.hypercube, symmetry)
#             throw(ArgumentError("lattice and point symmetry not compatible"))
#         end
#         elements = [embed(lattice, elem) for elem in symmetry.elements]
#         new(lattice, symmetry, elements)
#     end
# end


# elements(symbed::PointSymmetryEmbedding) = symbed.elements
# element(symbed::PointSymmetryEmbedding, g) = symbed.elements[g]


function little_symmetry(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                         psymbed::SymmetryEmbedding{PointSymmetry})
    tsymbed.lattice != psymbed.lattice && throw(ArgumentError("lattices do not match"))
    psym_little = little_symmetry(tsymbed.symmetry, psymbed.symmetry)
    return SymmetryEmbedding(psymbed.lattice, psym_little)
end

function little_symmetry(tsymbed::SymmetryEmbedding{TranslationSymmetry},
                         tsym_irrep::Integer,
                         psymbed::SymmetryEmbedding{PointSymmetry})
    psym_little = little_symmetry(tsymbed.symmetry, tsym_irrep, psymbed.symmetry)
    return SymmetryEmbedding(psymbed.lattice, psym_little)
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
