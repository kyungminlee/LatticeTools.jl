
export SymmorphicSymmetryEmbedding
export embed
export get_irrep_components

struct SymmorphicSymmetryEmbedding{S1<:AbstractSymmetry,
                                   S2<:AbstractSymmetry}<:AbstractSymmetryEmbedding
    lattice::Lattice
    normal::SymmetryEmbedding{S1}
    rest::SymmetryEmbedding{S2}

    function SymmorphicSymmetryEmbedding(lattice::Lattice,
            symmetry::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E}
        if !iscompatible(lattice, symmetry)
            throw(ArgumentError("lattice and symmetry are not compatible"))
        end
        return new{S1, S2}(lattice,
                           embed(lattice, symmetry.normal),
                           embed(lattice, symmetry.rest))
    end

    function SymmorphicSymmetryEmbedding(normal::SymmetryEmbedding{S1}, rest::SymmetryEmbedding{S2}) where {S1, S2}
        lattice = normal.lattice
        if rest.lattice != lattice
            throw(ArgumentError("two symmetry embeddings should have the same lattice"))
        end
        symmetry = SymmorphicSymmetry(normal.symmetry, rest.symmetry)
        return SymmorphicSymmetryEmbedding(lattice, symmetry)
    end
end

embed(lattice::Lattice, ssym::SymmorphicSymmetry) = SymmorphicSymmetryEmbedding(lattice, ssym)

function ⋊(normal::SymmetryEmbedding{S1}, rest::SymmetryEmbedding{S2}) where {S1, S2}
    lattice = normal.lattice
    if rest.lattice != lattice
        throw(ArgumentError("two symmetry embeddings should have the same lattice"))
    end
    return SymmorphicSymmetryEmbedding(lattice, SymmorphicSymmetry(normal.symmetry, rest.symmetry))
end

import Base.eltype
eltype(::SymmorphicSymmetryEmbedding) = SitePermutation
import Base.valtype
valtype(::SymmorphicSymmetryEmbedding) = SitePermutation

export iscompatible
function iscompatible(lattice::Lattice, ssym::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E}
    return iscompatible(lattice, ssym.normal) && iscompatible(lattice, ssym.rest)
end

# function get_irrep_components(sym::SymmorphicSymmetryEmbedding{S1, S2}) where {S1, S2}
#     (SymmorphicIrrepComponent(normal_sic.symmetry, normal_sic.irrep_index, normal_sic.irrep_component,
#                                    rest_sic.symmetry, rest_sic.irrep_index, rest_sic.irrep_component)
#         for normal_sic in get_irrep_components(sym.normal)
#         for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest)))
# end



