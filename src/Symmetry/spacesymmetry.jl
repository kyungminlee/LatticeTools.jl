

struct TranslationSymmetryEmbedding
    lattice::Lattice
    symmetry::TranslationSymmetry
    function TranslationSymmetryEmbedding(lattice::Lattice, symmetry::TranslationSymmetry)
        new(lattice, symmetry)
    end
end

struct PointSymmetryEmbedding
    lattice::Lattice
    symmetry::PointSymmetry

    function PointSymmetryEmbedding(lattice::Lattice, symmetry::PointSymmetry)
        if !iscompatible(lattice.hypercube, symmetry)
            throw(ArgumentError("lattice and point symmetry not compatible"))
        end
        new(lattice, symmetry)
    end
end

struct SymmorphicSpaceSymmetryEmbedding
    lattice::Lattice
    translation_symmetry::TranslationSymmetry
    point_symmetry::PointSymmetry

    function SymmorphicSpaceSymmetryEmbedding(
                lattice::Lattice,
                translation::TranslationSymmetry,
                point::PointSymmetry)
        if !iscompatible(translation, point)
            throw(ArgumentError("translation symmetry and point symmetry not compatible"))
        end
        new(lattice, translation, point)
    end
end
