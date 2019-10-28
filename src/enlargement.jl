
export make_supercell
export hypercubic_cluster
export ExactLinearAlgebra

using LinearAlgebra

module ExactLinearAlgebra

  function get_cofactor_matrix_unsafe!(out ::AbstractMatrix{I}, mat ::AbstractMatrix{I}, row ::Integer, col ::Integer) where {I<:Number}
    out[1:row-1, 1:col-1] = mat[1:row-1, 1:col-1]
    out[1:row-1, col:end] = mat[1:row-1, col+1:end]
    out[row:end, 1:col-1] = mat[row+1:end, 1:col-1]
    out[row:end, col:end] = mat[row+1:end, col+1:end]
    out
  end

  function determinant_unsafe(mat ::AbstractMatrix{S}) where {S<:Union{Integer,Complex{<:Integer}}}
    n = size(mat)[1]
    n == 1 && return mat[1]
    out = Matrix{S}(undef, (n-1, n-1))
    sign = 1
    D = zero(S)
    for f = 1:n
      get_cofactor_matrix_unsafe!(out, mat, 1, f)
      @inbounds D += sign * mat[1,f] * determinant_unsafe(out)
      sign = -sign
    end
    return D
  end

  function determinant(mat::AbstractMatrix{S}) where {S<:Union{Integer,Complex{<:Integer}}}
    n, m = size(mat)
    n != m && throw(ArgumentError("matrix needs to be square"))
    n <= 0 && throw(ArgumentError("matrix empty"))
    return determinant_unsafe(mat)
  end

  function inverse(mat::AbstractMatrix{S}) where {S<:Union{Integer,Complex{<:Integer}}}
    n, m = size(mat)
    n != m && throw(ArgumentError("matrix needs to be square"))
    n <= 0 && throw(ArgumentError("matrix empty"))
    n == 1 && return ones(S, (1,1)) .// mat[1]
    cofactor = Array{S}(undef, (n, n))
    temp = Array{S}(undef, (n-1, n-1))
    for r in 1:n, c in 1:n
      sign = 1 - 2 * mod(r+c, 2)
      get_cofactor_matrix_unsafe!(temp, mat, r, c)
      @inbounds cofactor[r,c] = sign * determinant_unsafe(temp)
    end
    D = sum(mat[1,:] .* cofactor[1,:])
    return transpose(cofactor) // D
  end

end # module ExactLinearAlgebra


export translation_element
function translation_element(displacement ::AbstractVector{<:Integer},
                             coordinates ::AbstractVector{<:AbstractVector{<:Integer}},
                             torus_wrap ::Function
                             )
  return Permutation([torus_wrap(r .+ displacement)[1] for r in coordinates])
end

export translation_symmetry_group
function translation_symmetry_group(hypercube ::HypercubicLattice)
  generator_indices = get_generators(hypercube)
  generators = [translation_element(hypercube, hypercube.coordinates[i]) for i in generator_indices]
  return TranslationGroup(generators)
end


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
