

struct TranslationSymmetryEmbedding
    lattice::Lattice
    symmetry::TranslationSymmetry
    operations::Vector{Permutation}

    function TranslationSymmetryEmbedding(lattice::Lattice, symmetry::TranslationSymmetry)
        permutations = get_orbital_permutations(lattice, symmetry)
        new(lattice, symmetry, permutations)
    end
end

struct PointSymmetryEmbedding
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

struct SymmorphicSpaceSymmetryEmbedding
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
