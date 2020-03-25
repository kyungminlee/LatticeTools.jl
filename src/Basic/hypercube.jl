export HypercubicLattice

export dimension
export isequiv
export translation_group_multiplication_table


struct HypercubicLattice
  scale_matrix ::Matrix{Int}

  inverse_scale_matrix ::Matrix{Rational{Int}}
  coordinates ::Vector{Vector{Int}}
  coordinate_indices ::Dict{Vector{Int}, Int}
  torus_wrap ::Function

  function HypercubicLattice(scale_matrix ::AbstractMatrix{<:Integer})
    n, m = size(scale_matrix)
    n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
    d = ExactLinearAlgebra.determinant(scale_matrix)
    d == 0 && throw(SingularException("scale matrix null"))
    d = abs(d)

    inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)

    max_range = sum((abs.(scale_matrix[:, i]) for i in 1:n))

    coords = Vector{Int}[]
    sizehint!(coords, abs(d))
    for g in Iterators.product([0:2*x+1 for x in max_range]...)
      r1 = collect(g)
      r2 = inverse_scale_matrix * r1
      R = Int.(floor.(r2))
      r3 = r1 - scale_matrix * R
      if !(r3 in coords)
        push!(coords, r3)
      end
    end
    @assert length(coords) == d
    coord_indices = Dict(r => i for (i, r) in enumerate(coords))

    function torus_wrap(r ::AbstractVector{<:Integer})
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, coord_indices[r2]
    end

    function torus_wrap(r ::AbstractMatrix{<:Integer})
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R # TODO: need to be tested
      return R, [coord_indices[x] for x in eachcol(r2)]
    end

    @assert all(torus_wrap(r) == (zeros(n), i) for (i, r) in enumerate(coords))

    return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, torus_wrap)
  end
end


dimension(hypercube::HypercubicLattice) = size(hypercube.scale_matrix, 1)


function isequiv(lhs::HypercubicLattice, rhs::HypercubicLattice)
    det_lhs = TightBindingLattice.ExactLinearAlgebra.determinant(lhs.scale_matrix)
    det_rhs = TightBindingLattice.ExactLinearAlgebra.determinant(rhs.scale_matrix)
    det_lhs != det_rhs && return false

    inv_lhs = TightBindingLattice.ExactLinearAlgebra.inverse(lhs.scale_matrix)
    inv_rhs = TightBindingLattice.ExactLinearAlgebra.inverse(rhs.scale_matrix)

    return all(isinteger.(inv_lhs * rhs.scale_matrix)) && all(isinteger.(inv_rhs * lhs.scale_matrix))
end

function translation_group_multiplication_table(hypercube::HypercubicLattice)
    n = length(hypercube.coordinates)
    mtab = zeros(Int, (n, n))
    wrap = hypercube.torus_wrap
    for (i, ri) in enumerate(hypercube.coordinates), (j, rj) in enumerate(hypercube.coordinates)
        _, k = wrap(ri.+rj)
        mtab[i,j] = k
    end
    return mtab
end
