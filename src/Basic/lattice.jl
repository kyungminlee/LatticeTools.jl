export make_lattice, makelattice
export Lattice
export dimension


"""
    Lattice{O}

Represent a lattice.

# Arguments
* `unitcell::UnitCell{O}`
* `hypercube::Hypercube`
* `bravais_coordinates::Vector{Vector{Int}}`
* `supercell::UnitCell{Tuple{O, Vector{Int}}}`
"""
struct Lattice{O}
    unitcell::UnitCell{O}
    hypercube::Hypercube
    bravais_coordinates::Vector{Vector{Int}}
    supercell::UnitCell{Tuple{O, Vector{Int}}}
end


function Base.:(==)(lhs::Lattice{O}, rhs::Lattice{O}) where O
    return (
        lhs.unitcell == rhs.unitcell
        && lhs.hypercube == rhs.hypercube
        && lhs.bravais_coordinates == rhs.bravais_coordinates
    )
end


function makelattice(
    unitcell::UnitCell{O},
    shape_matrix::AbstractMatrix{<:Integer}
) where O
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end

    hypercube = Hypercube(shape_matrix)
    generator_translations = find_generators(hypercube)
    unitcell_coordinates = generate_coordinates(hypercube, generator_translations)

    new_latticevectors = unitcell.latticevectors * hypercube.shape_matrix

    new_unitcell = makeunitcell(new_latticevectors; SiteType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.sites
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addsite!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, hypercube, unitcell_coordinates, new_unitcell)
end
#TODO unit testing for lattice with wrong dimensions


"""
    makelattice(unitcell, shape, [generators])

Create a lattice with periodic boundary condition, using the unitcell, shape, and
translation generators.

# Arguments
* `unitcell::UnitCell{O}`
* `shape_matrix::AbstractMatrix{<:Integer}`: shape of the Bravais lattice. This can also be
a single integer, which is equivalent to a identity matrix times the number.

# Optional Arguments
* `generator_translations::AbstractMatrix{<:Integer}`
"""
function makelattice(
    unitcell::UnitCell{O},
    shape_matrix::AbstractMatrix{<:Integer},
    generator_translations::AbstractMatrix{<:Integer},
) where O
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end
    hypercube = Hypercube(shape_matrix)
    unitcell_coordinates = generate_coordinates(hypercube, generator_translations)

    new_latticevectors = unitcell.latticevectors * hypercube.shape_matrix

    new_unitcell = makeunitcell(new_latticevectors; SiteType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.sites
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addsite!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, hypercube, unitcell_coordinates, new_unitcell)
end


function makelattice(unitcell::UnitCell, scale::Integer)
    return makelattice(unitcell, scale*ones(Int, (1,1)))
end


make_lattice = makelattice

"""
    dimension(lattice)

Spatial dimension of the lattice.

# Arguments
- `lattice::Lattice`
"""
dimension(lattice::Lattice) = dimension(lattice.unitcell)
