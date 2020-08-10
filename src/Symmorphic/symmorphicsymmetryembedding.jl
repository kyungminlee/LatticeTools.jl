
export SymmorphicSymmetryEmbedding
export embed
export get_irrep_components
export iscompatible
export ⋉, ⋊

export fractional_momentum

struct SymmorphicSymmetryEmbedding{S1<:AbstractSymmetry,
                                   S2<:AbstractSymmetry}<:AbstractSymmetryEmbedding
    lattice::Lattice
    normal::SymmetryEmbedding{S1}
    rest::SymmetryEmbedding{S2}

    function SymmorphicSymmetryEmbedding(
        lattice::Lattice,
        symmetry::SymmorphicSymmetry{S1, S2, E},
    ) where {S1, S2, E}
        if !iscompatible(lattice, symmetry)
            throw(ArgumentError("lattice and symmetry are not compatible"))
        end
        return new{S1, S2}(
            lattice,
            embed(lattice, symmetry.normal),
            embed(lattice, symmetry.rest),
        )
    end

    function SymmorphicSymmetryEmbedding(normal::SymmetryEmbedding, rest::SymmetryEmbedding)
        lattice = normal.lattice
        if rest.lattice != lattice
            throw(ArgumentError("two symmetry embeddings should have the same lattice"))
        end
        symmetry = SymmorphicSymmetry(normal.symmetry, rest.symmetry)
        return SymmorphicSymmetryEmbedding(lattice, symmetry)
    end
end


"""
    embed(lattice::Lattice, ssym::SymmorphicSymmetry)
"""
function embed(lattice::Lattice, ssym::SymmorphicSymmetry)
    return SymmorphicSymmetryEmbedding(lattice, ssym)
end


function ⋊(normal::SymmetryEmbedding, rest::SymmetryEmbedding)
    return SymmorphicSymmetryEmbedding(normal, rest)
end

function ⋉(rest::SymmetryEmbedding, normal::SymmetryEmbedding)
    return SymmorphicSymmetryEmbedding(normal, rest)
end


Base.eltype(::SymmorphicSymmetryEmbedding) = SitePermutation

Base.valtype(::SymmorphicSymmetryEmbedding) = SitePermutation


function Base.:(==)(lhs::SSE, rhs::SSE) where {SSE<:SymmorphicSymmetryEmbedding}
    return lhs.lattice == rhs.lattice && lhs.normal == rhs.normal && lhs.rest == rhs.rest
end


group_order(sym::SymmorphicSymmetryEmbedding) = group_order(sym.normal) * group_order(sym.rest)


function fractional_momentum(sym::SymmorphicSymmetryEmbedding{<:TranslationSymmetry, S2}, args...) where S2
    return fractional_momentum(sym.normal, args...)
end


"""
    iscompatible(lattice::Lattice, ssym::SymmorphicSymmetry)
"""
function iscompatible(lattice::Lattice, ssym::SymmorphicSymmetry)
    return iscompatible(lattice, ssym.normal) && iscompatible(lattice, ssym.rest)
end
