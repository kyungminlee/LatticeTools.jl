export AbstractSymmetryEmbedding

export TranslationSymmetryEmbedding
export PointSymmetryEmbedding
export SymmorphicSpaceSymmetryEmbedding

export embed

abstract type AbstractSymmetryEmbedding end
abstract type AbstractSymmetryOperationEmbedding <:AbstractSymmetryOperation end


struct ProductSymmetry{S1<:AbstractSymmetry, S2<:AbstractSymmetry}
    symmetry1::S1
    symmetry2::S2

    function ProductSymmetry(sym1::S1, sym2::S2) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        new{S1, S2}(sym1, sym2)
    end
end


struct TranslationSymmetryEmbedding<:AbstractSymmetryEmbedding
    lattice::Lattice
    symmetry::TranslationSymmetry
    operations::Vector{Permutation}

    function TranslationSymmetryEmbedding(lattice::Lattice, symmetry::TranslationSymmetry)
        if (lattice.hypercube != symmetry.hypercube)
            throw(ArgumentError("lattice and translation symmetry not compatible"))
        end
        permutations = get_orbital_permutations(lattice, symmetry)
        new(lattice, symmetry, permutations)
    end
end

struct PointSymmetryEmbedding<:AbstractSymmetryEmbedding
    lattice::Lattice
    symmetry::PointSymmetry
    operations::Vector{Permutation}

    function PointSymmetryEmbedding(lattice::Lattice, symmetry::PointSymmetry)
        if !iscompatible(lattice.hypercube, symmetry)
            throw(ArgumentError("lattice and point symmetry not compatible"))
        end
        permutations = get_orbital_permutations(lattice, symmetry)
        new(lattice, symmetry, permutations)
    end
end


struct SymmorphicSpaceSymmetryEmbedding<:AbstractSymmetryEmbedding
    lattice::Lattice
    translation_symmetry::TranslationSymmetry
    point_symmetry::PointSymmetry
    operations::Matrix{Permutation}

    function SymmorphicSpaceSymmetryEmbedding(
                lattice::Lattice,
                translation::TranslationSymmetry,
                point::PointSymmetry)
        if !iscompatible(translation, point)
            throw(ArgumentError("translation symmetry and point symmetry not compatible"))
        end
        tperms = get_orbital_permutations(lattice, translation)
        pperms = get_orbital_permutations(lattice, point)
        permutations = [pop * top for top in tperms, pop in pperms]
        @assert size(permutations) == (length(tperms), length(pperms))
        new(lattice, translation, point, permutations)
    end
end



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
        ($f)(symbed::TranslationSymmetryEmbedding, args...) = ($f)(symbed.symmetry, args...)
        ($f)(symbed::PointSymmetryEmbedding, args...) = ($f)(symbed.symmetry, args...)
    end)
end


function embed(lattice::Lattice, tsym::TranslationSymmetry)
    TranslationSymmetryEmbedding(lattice, tsym)
end

function embed(lattice::Lattice, psym::PointSymmetry)
    PointSymmetryEmbedding(lattice, psym)
end

function embed(lattice::Lattice, tsym::TranslationSymmetry, psym::PointSymmetry)
    SymmorphicSpaceSymmetryEmbedding(lattice, tsym, psym)
end
