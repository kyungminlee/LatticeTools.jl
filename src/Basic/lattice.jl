export make_lattice, makelattice
export Lattice
export dimension


"""
    Lattice{S, O}

Represent a lattice.

# Arguments
* `unitcell::UnitCell{S, O}`
* `hypercube::Hypercube`
* `bravais_coordinates::Vector{Vector{Int}}`
* `supercell::UnitCell{Tuple{S, Vector{Int}}, Tuple{O, Vector{Int}}}`
"""
struct Lattice{S, O}
    unitcell::UnitCell{S, O}
    hypercube::Hypercube
    bravais_coordinates::Vector{Vector{Int}}
    supercell::UnitCell{Tuple{S, Vector{Int}}, Tuple{O, Vector{Int}}}
end


function Base.:(==)(lhs::Lattice{S, O}, rhs::Lattice{S, O}) where {S, O}
    return (
        lhs.unitcell == rhs.unitcell
        && lhs.hypercube == rhs.hypercube
        && lhs.bravais_coordinates == rhs.bravais_coordinates
    )
end


"""
    makelattice(unitcell, shape)

Create a `Lattice` from `unitcell` and `shape`.
Generating translations are fouond using `find_generators`.

# Arguments
- `unitcell`
- `shape`: shape of the cluster
"""
function makelattice(
    unitcell::UnitCell{S, O},
    shape_matrix::AbstractMatrix{<:Integer}
) where {S, O}
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end

    hypercube = Hypercube(shape_matrix)
    generator_translations = find_generators(hypercube)
    unitcell_coordinates = generate_coordinates(hypercube, generator_translations)

    new_latticevectors = unitcell.latticevectors * hypercube.shape_matrix

    new_unitcell = makeunitcell(new_latticevectors; SiteType=Tuple{S, Vector{Int}}, OrbitalType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates
        for (sitename, sitecoord) in unitcell.sites
            cc = fract2carte(unitcell, sitecoord)
            new_cc = cc + unitcell.latticevectors * uc
            new_sitecoord = carte2fract(new_unitcell, new_cc)
            new_sitename = (sitename, uc)
            addsite!(new_unitcell, new_sitename, new_sitecoord)
        end
        for (orbitalname, orbitalsiteindex) in unitcell.orbitals
            sitename = getsitename(unitcell, orbitalsiteindex)
            new_orbitalname = (orbitalname, uc)
            new_sitename = (sitename, uc)
            new_orbitalsiteindex = getsiteindex(new_unitcell, new_sitename)
            addorbital!(new_unitcell, new_orbitalname, new_orbitalsiteindex)
        end
    end
    return Lattice{S, O}(unitcell, hypercube, unitcell_coordinates, new_unitcell)
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
    unitcell::UnitCell{S, O},
    shape_matrix::AbstractMatrix{<:Integer},
    generator_translations::AbstractMatrix{<:Integer},
) where {S, O}
    dim = dimension(unitcell)
    if size(shape_matrix) != (dim, dim)
        throw(DimensionMismatch("unitcell and shape_matrix should have the same dimension"))
    end
    hypercube = Hypercube(shape_matrix)
    unitcell_coordinates = generate_coordinates(hypercube, generator_translations)

    new_latticevectors = unitcell.latticevectors * hypercube.shape_matrix

    new_unitcell = makeunitcell(new_latticevectors; SiteType=Tuple{S, Vector{Int}}, OrbitalType=Tuple{O, Vector{Int}})
    for uc in unitcell_coordinates
        for (sitename, sitecoord) in unitcell.sites
            cc = fract2carte(unitcell, sitecoord)
            new_cc = cc + unitcell.latticevectors * uc
            new_sitecoord = carte2fract(new_unitcell, new_cc)
            new_sitename = (sitename, uc)
            addsite!(new_unitcell, new_sitename, new_sitecoord)
        end
        for (orbitalname, orbitalsiteindex) in unitcell.orbitals
            sitename = getsitename(unitcell, orbitalsiteindex)
            new_orbitalname = (orbitalname, uc)
            new_sitename = (sitename, uc)
            new_orbitalsiteindex = getsiteindex(new_unitcell, new_sitename)
            addorbital!(new_unitcell, new_orbitalname, new_orbitalsiteindex)
        end
    end

    return Lattice{S, O}(unitcell, hypercube, unitcell_coordinates, new_unitcell)
end


"""
    makelattice(unitcell, scale)

Construct a `Lattice` with `unitcell` having shape `scale`.
Works for one-dimensional lattice.
"""
function makelattice(unitcell::UnitCell, scale::Integer)
    return makelattice(unitcell, hcat(scale), hcat(1))
end


const make_lattice = makelattice


"""
    dimension(lattice)

Spatial dimension of the lattice.

# Arguments
- `lattice::Lattice`
"""
dimension(lattice::Lattice) = dimension(lattice.unitcell)
