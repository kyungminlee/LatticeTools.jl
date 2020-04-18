export make_lattice
export Lattice
export dimension


struct Lattice{O}
    unitcell::UnitCell{O}
    hypercube::HypercubicLattice
    supercell::UnitCell{Tuple{O, Vector{Int}}}
end

import Base.==
function ==(lhs::Lattice{O}, rhs::Lattice{O}) where O
    return lhs.unitcell == rhs.unitcell && lhs.hypercube == rhs.hypercube
end


function make_lattice(unitcell::UnitCell, scale::Integer)
    make_lattice(unitcell, scale*ones(Int, (1,1)))
end

function make_lattice(unitcell::UnitCell{O}, scale_matrix::AbstractMatrix{<:Integer}) where O
    dim = dimension(unitcell)
    if size(scale_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and scale_matrix should have the same dimension"))
    end
    
    hypercube = orthogonalize(HypercubicLattice(scale_matrix))

    new_latticevectors = unitcell.latticevectors * hypercube.scale_matrix
    inverse_scale_matrix = hypercube.inverse_scale_matrix
    unitcell_coordinates = hypercube.coordinates

    new_unitcell = make_unitcell(new_latticevectors; OrbitalType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.orbitals
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addorbital!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, hypercube, new_unitcell)
end
