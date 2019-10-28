
export make_supercell
export hypercubic_cluster



#=
 . . . . . o . .
 . . . o . . . .
 . o . . . . o .
 . . . . o . . .
 . . o V . . . o
 o . . . . o . .
 . . . o . . . .
 . o . . . . . .
=#
"""

# Returns
  - new unit cell
  - embedding function
    - which takes the unitcell displacement in the orignal displacement
    - which returns a tuple of (supercell lattice displacement, sublattice subscript)
"""
function make_supercell(unitcell ::UnitCell{O}, scale_matrix ::AbstractMatrix{<:Integer}) where O
  # check dimensions
  new_latticevectors = unitcell.latticevectors * scale_matrix
  inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)
  unitcell_coordinates, torus_wrap = hypercubic_cluster(scale_matrix; inverse_scale_matrix=inverse_scale_matrix)
  new_unitcell = make_unitcell(new_latticevectors; OrbitalType=Tuple{O, Vector{Int}})
  for uc in unitcell_coordinates, (orbname, orbcoord) in unitcell.orbitals
    cc = fract2carte(unitcell, orbcoord)
    new_cc = cc + unitcell.latticevectors * uc
    new_orbcoord = carte2fract(new_unitcell, new_cc)
    new_orbname = (orbname, uc)
    addorbital!(new_unitcell, new_orbname, new_orbcoord)
  end

  # function torus_wrap(lattice_displacement ::AbstractVector{<:Integer})
  #   R2 = I.(floor.(inverse_scale_matrix * lattice_displacement))
  #   r2 = lattice_displacement - scale_matrix * R2
  #   return R2, r2
  # end

  return (new_unitcell, torus_wrap)
end
