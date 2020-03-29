export HypercubicLattice

export dimension
export isequiv
export translation_group_multiplication_table
export orthogonalize


struct HypercubicLattice
  scale_matrix ::Matrix{Int}

  inverse_scale_matrix ::Matrix{Rational{Int}}
  coordinates ::Vector{Vector{Int}}
  coordinate_indices ::Dict{Vector{Int}, Int}
  torus_wrap ::Function
  wrap ::Function

  function HypercubicLattice(scale_matrix ::AbstractMatrix{<:Integer})
    n, m = size(scale_matrix)
    n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
    d = ExactLinearAlgebra.determinant(scale_matrix)
    d == 0 && throw(SingularException("scale matrix null"))
    d = abs(d)

    inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)

    function wrap(r ::AbstractVector{<:Integer})
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, r2
    end

    function wrap(r ::AbstractMatrix{<:Integer})
      @warn "Need to be tested"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, r2
    end

    max_range = sum((abs.(scale_matrix[:, i]) for i in 1:n))
    coords = Vector{Int}[]
    sizehint!(coords, abs(d))
    for g in Iterators.product([0:2*x+1 for x in max_range]...)
      r1 = collect(g)
      _, r2 = wrap(r1)
      if !(r2 in coords)
        push!(coords, r2)
      end
    end

    @assert length(coords) == d
    coord_indices = Dict(r => i for (i, r) in enumerate(coords))

    function torus_wrap(r ::AbstractVector{<:Integer})
      @warn "torus_wrap is deprecated"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, coord_indices[r2]
    end

    function torus_wrap(r ::AbstractMatrix{<:Integer})
      @warn "torus_wrap is deprecated"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R # TODO: need to be tested
      return R, [coord_indices[x] for x in eachcol(r2)]
    end

    @assert all(wrap(r) == (zeros(n), r) for (i, r) in enumerate(coords))

    return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, torus_wrap, wrap)
  end


  function HypercubicLattice(scale_matrix ::AbstractMatrix{<:Integer},
                             coords::AbstractVector{<:AbstractVector{<:Integer}})
    n, m = size(scale_matrix)
    n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
    d = ExactLinearAlgebra.determinant(scale_matrix)
    d == 0 && throw(SingularException("scale matrix null"))
    d = abs(d)

    inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)

    coord_indices = Dict{Vector{Int}, Int}(r => i for (i, r) in enumerate(coords))

    function torus_wrap(r ::AbstractVector{<:Integer})
      @warn "torus_wrap is deprecated"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, coord_indices[r2]
    end

    function torus_wrap(r ::AbstractMatrix{<:Integer})
      @warn "torus_wrap is deprecated"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R # TODO: need to be tested
      return R, [coord_indices[x] for x in eachcol(r2)]
    end

    function wrap(r ::AbstractVector{<:Integer})
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return R, r2
    end

    function wrap(r ::AbstractMatrix{<:Integer})
      @warn "Need to be tested"
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R # TODO: need to be tested
      return R, r2
    end

    if !all(wrap(r) == (zeros(n), r) for (i, r) in enumerate(coords))
      throw(ArgumentError("coordinates are not that of the hypercubic lattice with size $(scale_matrix)"))
    end
    !allunique(coords) && throw(ArgumentError("coordinates not unique"))
    length(coords) != d && throw(ArgumentError("too few coordinates"))

    return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, torus_wrap, wrap)
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
    wrap = hypercube.wrap
    for (i, ri) in enumerate(hypercube.coordinates), (j, rj) in enumerate(hypercube.coordinates)
        _, rk = wrap(ri.+rj)
        k = hypercube.coordinate_indices[rk]
        mtab[i,j] = k
    end
    return mtab
end
