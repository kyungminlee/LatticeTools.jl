export make_lattice, makelattice
export Lattice
export dimension


"""
    Lattice{O}

Represent a lattice.

# Arguments
* `unitcell::UnitCell{O}`
* `orthocube::OrthoCube`
* `bravais_coordinates::Vector{Vector{Int}}`
* `supercell::UnitCell{Tuple{O, Vector{Int}}}`
"""
struct Lattice{O}
    unitcell::UnitCell{O}
    orthocube::OrthoCube
    bravais_coordinates::Vector{Vector{Int}}
    supercell::UnitCell{Tuple{O, Vector{Int}}}
end


function Base.:(==)(lhs::Lattice{O}, rhs::Lattice{O}) where O
    return (
        lhs.unitcell == rhs.unitcell
        && lhs.orthocube == rhs.orthocube
        && lhs.bravais_coordinates == rhs.bravais_coordinates
    )
end


function make_lattice(unitcell::UnitCell, scale::Integer)
    return make_lattice(unitcell, scale*ones(Int, (1,1)))
end

function make_lattice(
    unitcell::UnitCell{O},
    shape_matrix::AbstractMatrix{<:Integer}
) where O
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end

    orthocube = OrthoCube(shape_matrix)
    generator_translations = find_generators(orthocube)
    unitcell_coordinates = generate_coordinates(orthocube, generator_translations)

    new_latticevectors = unitcell.latticevectors * orthocube.shape_matrix

    new_unitcell = make_unitcell(new_latticevectors; SiteType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.sites
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addsite!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, orthocube, unitcell_coordinates, new_unitcell)
end
#TODO unit testing for lattice with wrong dimensions


function make_lattice(
    unitcell::UnitCell{O},
    shape_matrix::AbstractMatrix{<:Integer},
    generator_translations::AbstractMatrix{<:Integer},
) where O
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end
    orthocube = OrthoCube(shape_matrix)
    unitcell_coordinates = generate_coordinates(orthocube, generator_translations)

    new_latticevectors = unitcell.latticevectors * orthocube.shape_matrix

    new_unitcell = make_unitcell(new_latticevectors; SiteType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.sites
        cc = fract2carte(unitcell, orbcoord)
        new_cc = cc + unitcell.latticevectors * uc
        new_orbcoord = carte2fract(new_unitcell, new_cc)
        new_orbname = (orbname, uc)
        addsite!(new_unitcell, new_orbname, new_orbcoord)
    end

    return Lattice{O}(unitcell, orthocube, unitcell_coordinates, new_unitcell)
end

makelattice = make_lattice

dimension(lattice::Lattice) = dimension(lattice.unitcell)
