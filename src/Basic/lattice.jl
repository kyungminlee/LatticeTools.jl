export SuperCell
export make_supercell
export make_lattice
export Lattice


#=
struct SuperCell{O}
    unitcell::UnitCell{O}
    hypercube::HypercubicLattice
    supercell::UnitCell{Tuple{O, Vector{Int}}} # The supercell as a single unit cell
end


function make_supercell(unitcell ::UnitCell{O}, scale_matrix ::AbstractMatrix{<:Integer}) where O
    hypercube = HypercubicLattice(scale_matrix)
    return make_supercell(unitcell, hypercube)
end


function make_supercell(unitcell ::UnitCell{O}, hypercube::HypercubicLattice) where O
  # check dimensions
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

  return SuperCell{O}(unitcell, hypercube, new_unitcell)
end
=#


struct Lattice{O}
  unitcell::UnitCell{O}
  hypercube::HypercubicLattice
  supercell::UnitCell{Tuple{O, Vector{Int}}}
end

function make_lattice(unitcell::UnitCell, scale::Integer)
  make_lattice(unitcell, scale*ones(Int, (1,1)))
end


function make_lattice(unitcell::UnitCell{O}, scale_matrix::AbstractMatrix{<:Integer}) where O
  hypercube = orthogonalize(HypercubicLattice(scale_matrix))
  #return make_lattice(unitcell, hypercube)
#end

#function make_lattice(unitcell::UnitCell{O}, hypercube::HypercubicLattice) where O
  #hypercube = orthogonalize(hypercube)

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



#=
save for now
"""
    make_supercell(unitcell, scale_matrix; compute_symmetry=true)

# Arguments
  - `unitcell`
  - `scale_matrix`
  - `compute_symmetry=true` (optional)
# Returns
  - new unit cell
  - embedding function
    - which takes the unitcell displacement in the orignal displacement
    - which returns a tuple of (supercell lattice displacement, sublattice subscript)
  - translation group
"""
function make_supercell(unitcell ::UnitCell{O},
                        scale_matrix ::AbstractMatrix{<:Integer};
                        compute_symmetry::Bool=true) where O
  # check dimensions
  new_latticevectors = unitcell.latticevectors * scale_matrix
  hypercube = HypercubicLattice(scale_matrix)
  inverse_scale_matrix = hypercube.inverse_scale_matrix
  unitcell_coordinates = hypercube.coordinates
  torus_wrap = hypercube.torus_wrap

  new_unitcell = make_unitcell(new_latticevectors; OrbitalType=Tuple{O, Vector{Int}})
  for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.orbitals
    cc = fract2carte(unitcell, orbcoord)
    new_cc = cc + unitcell.latticevectors * uc
    new_orbcoord = carte2fract(new_unitcell, new_cc)
    new_orbname = (orbname, uc)
    addorbital!(new_unitcell, new_orbname, new_orbcoord)
  end

  if compute_symmetry
    translation_group = translation_symmetry_group(hypercube)
    new_generators = Permutation[]
    for p in translation_group.generators
        p2 = Vector{Int}(undef, length(p.map) * numorbital(unitcell))
        for (i, j) in enumerate(p.map)
            ri = hypercube.coordinates[i]
            rj = hypercube.coordinates[j]
            for (orbname, orbcoord) in unitcell.orbitals
                ip = getorbitalindex(new_unitcell, (orbname, ri))
                jp = getorbitalindex(new_unitcell, (orbname, rj))
                p2[ip] = jp
            end
        end
        push!(new_generators, Permutation(p2))
    end

    return (new_unitcell, hypercube, TranslationGroup(new_generators))
  else
    return (new_unitcell, hypercube, nothing)
  end
end
=#
