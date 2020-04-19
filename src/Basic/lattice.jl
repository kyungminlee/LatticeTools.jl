export make_lattice
export Lattice
export dimension


struct Lattice{O}
    unitcell::UnitCell{O}
    orthocube::OrthoCube
    bravais_coordinates::Vector{Vector{Int}}
    supercell::UnitCell{Tuple{O, Vector{Int}}}
end

import Base.==
function ==(lhs::Lattice{O}, rhs::Lattice{O}) where O
    return lhs.unitcell == rhs.unitcell && lhs.orthocube == rhs.orthocube
end


function make_lattice(unitcell::UnitCell, scale::Integer)
    make_lattice(unitcell, scale*ones(Int, (1,1)))
end

function make_lattice(unitcell::UnitCell{O}, scale_matrix::AbstractMatrix{<:Integer}) where O
    dim = dimension(unitcell)
    if size(scale_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and scale_matrix should have the same dimension"))
    end
    
    # hypercube = orthogonalize(HypercubicLattice(scale_matrix))
    orthocube = OrthoCube(scale_matrix)
    generator_translations = find_generators(orthocube)
    coordinates = generate_coordinates(orthocube, generator_translations)

    # new_latticevectors = unitcell.latticevectors * hypercube.scale_matrix
    # inverse_scale_matrix = hypercube.inverse_scale_matrix
    # unitcell_coordinates = hypercube.coordinates
    
    new_latticevectors = unitcell.latticevectors * orthocube.shape_matrix
    inverse_scale_matrix = orthocube.inverse_shape_matrix
    unitcell_coordinates = coordinates

    new_unitcell = make_unitcell(new_latticevectors; OrbitalType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.orbitals
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addorbital!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, orthocube, unitcell_coordinates, new_unitcell)
end
#TODO unit testing for lattice with wrong dimensions

dimension(lattice::Lattice) = dimension(lattice.unitcell)