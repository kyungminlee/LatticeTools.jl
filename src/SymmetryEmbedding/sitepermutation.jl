export SitePermutation
export embed
export inverse

struct SitePermutation <:AbstractSymmetryOperationEmbedding
    permutation::Permutation
    SitePermutation(p::Permutation) = new(p)
    SitePermutation(p::AbstractVector{<:Integer}) = new(Permutation(p))
end


import Base.==
==(lhs::SitePermutation, rhs::SitePermutation) = lhs.permutation == rhs.permutation

import Base.isequal
isequal(lhs::SitePermutation, rhs::SitePermutation) = isequal(lhs.permutation, rhs.permutation)

import Base.*
*(lhs::SitePermutation, rhs::SitePermutation) = SitePermutation(lhs.permutation * rhs.permutation)

import Base.^
^(lhs::SitePermutation, rhs::Integer) = SitePermutation(lhs.permutation^rhs)

import Base.isless
isless(lhs::SitePermutation, rhs::SitePermutation) = isless(lhs.permutation, rhs.permutation)

import Base.hash
hash(arg::SitePermutation) = hash(arg.permutation)


inverse(sp::SitePermutation) = SitePermutation(inverse(sp.permutation))



function embed(lattice::Lattice, top::TranslationOperation)
    p = zeros(Int, numorbital(lattice.supercell))
    for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
        _, uc_coord2 = lattice.hypercube.wrap( top(uc_coord1) )
        orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
        orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
        p[orbital_index1] = orbital_index2
    end
    return SitePermutation(p)
end


function embed(lattice::Lattice, pop::PointOperation,
               orbital_map::AbstractVector{<:Tuple{<:Integer, <:AbstractVector{<:Integer}}})
    p = zeros(Int, numorbital(lattice.supercell))
    for (i, (j, dR)) in enumerate(orbital_map)
        namei = getorbitalname(lattice.unitcell, i)
        namej = getorbitalname(lattice.unitcell, j)
        for Ri in lattice.hypercube.coordinates
            _, Rj = lattice.hypercube.wrap(pop.matrix * Ri + dR)
            i_super = lattice.supercell.orbitalindices[(namei, Ri)]
            j_super = lattice.supercell.orbitalindices[(namej, Rj)]
            p[i_super] = j_super
        end
    end
    return SitePermutation(p)
end

