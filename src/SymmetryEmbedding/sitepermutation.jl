export SitePermutation
export embed
export isidentity


abstract type AbstractSymmetryOperationEmbedding <: AbstractSymmetryOperation{Int} end

"""
    SitePermutation
"""
struct SitePermutation <:AbstractSymmetryOperationEmbedding
    permutation::Permutation
    SitePermutation(p::Permutation) = new(p)
    SitePermutation(p::AbstractVector{<:Integer}) = new(Permutation(p))
end


import Base.==
(==)(lhs::SitePermutation, rhs::SitePermutation) = lhs.permutation == rhs.permutation

import Base.*
(*)(lhs::SitePermutation, rhs::SitePermutation) = SitePermutation(lhs.permutation * rhs.permutation)

import Base.^
(^)(lhs::SitePermutation, rhs::Integer) = SitePermutation(lhs.permutation^rhs)

# import Base.isless
# isless(lhs::SitePermutation, rhs::SitePermutation) = isless(lhs.permutation, rhs.permutation)

import Base.hash

function hash(arg::SitePermutation, h::UInt=UInt(0x0))
    return hash(arg.permutation, hash(SitePermutation, h))
end

import Base.inv
inv(sp::SitePermutation) = SitePermutation(inv(sp.permutation))

(p::SitePermutation)(i::Integer) = p.permutation(i)


"""
    embed(lattice, translation_operation)

Embed the simplest version of integer translation (no mapping between sites etc.)
"""
function embed(lattice::Lattice, top::TranslationOperation{<:Integer})
    p = zeros(Int, numsite(lattice.supercell))
    for (site_index1, ((site_name1, uc_coord1), _)) in enumerate(lattice.supercell.sites)
        _, uc_coord2 = lattice.orthocube.wrap( top(uc_coord1) )
        site_index1 = getsiteindex(lattice.supercell, (site_name1, uc_coord1))
        site_index2 = getsiteindex(lattice.supercell, (site_name1, uc_coord2))
        p[site_index1] = site_index2
    end
    return SitePermutation(p)
end


"""
    embed(lattice, point_operation)

Embed the simplest version of point operation. (no local unitary operation)
"""
function embed(lattice::Lattice, pop::PointOperation{<:Integer})
    site_map = findsitemap(lattice.unitcell, pop)
    if isnothing(site_map)
        throw(ArgumentError("lattice not compatible with $pop"))
    end
    p = zeros(Int, numsite(lattice.supercell))
    for (i, (j, dR)) in enumerate(site_map)
        namei = getsitename(lattice.unitcell, i)
        namej = getsitename(lattice.unitcell, j)
        for Ri in lattice.bravais_coordinates
            _, Rj = lattice.orthocube.wrap(pop.matrix * Ri + dR)
            i_super = lattice.supercell.siteindices[(namei, Ri)]
            j_super = lattice.supercell.siteindices[(namej, Rj)]
            p[i_super] = j_super
        end
    end
    return SitePermutation(p)
end


"""
    embed(lattice, sop::SpaceOperation{<:Integer, <:Integer})
"""
function embed(lattice::Lattice, sop::SpaceOperation{<:Integer, <:Integer})
    embed(lattice, PointOperation(sop.matrix)) * embed(lattice, TranslationOperation(sop.displacement))
end


isidentity(perm::SitePermutation) = isidentity(perm.permutation)