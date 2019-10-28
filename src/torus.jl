export HypercubicLattice
export get_generators

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

    inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)

    max_range = sum((abs.(scale_matrix[:, i]) for i in 1:n))

    coords = Vector{Int}[]
    sizehint!(coords, d)
    for g in Iterators.product([0:2*x+1 for x in max_range]...)
      #for g in Iterators.product(reverse([-x:x for x in reverse(max_range)])...)
      r1 = collect(g)
      r2 = inverse_scale_matrix * r1
      R = Int.(floor.(r2))
      r3 = r1 - scale_matrix * R
      #if all(0 <= x < 1 for x in r2)
      if !(r3 in coords)
        push!(coords, r3)
      end
    end
    @assert length(coords) == d
    coord_indices = Dict(r => i for (i, r) in enumerate(coords))
    function torus_wrap(r ::AbstractVector{<:Integer})
      R = Int.(floor.(inverse_scale_matrix * r))
      r2 = r - scale_matrix * R
      return coord_indices[r2], R
    end
    @assert all(torus_wrap(r) == (i, zeros(n)) for (i, r) in enumerate(coords))

    return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, torus_wrap)
  end
end

"""
    hypercubic_cluster

Generate a hypercubic cluster

```
 . . . . . .
 . . . o . .
 . o . . . .
 . . . . o .
 . . o . . .
 . . . . . .
```
"""
function hypercubic_cluster(scale_matrix ::AbstractMatrix{<:Integer}; inverse_scale_matrix ::Union{Nothing, AbstractMatrix{<:Rational}, AbstractMatrix{<:Integer}}=nothing)
  n, m = size(scale_matrix)
  n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
  d = ExactLinearAlgebra.determinant(scale_matrix)
  d == 0 && throw(SingularException("scale matrix null"))

  if inverse_scale_matrix === nothing
    inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)
  else
    (inverse_scale_matrix * scale_matrix == LinearAlgebra.I) || throw(ArgumentError("inverse_scale_matrix is not the inverse of scale_matrix"))
  end
  max_range = sum((abs.(scale_matrix[:, i]) for i in 1:n))

  coords = Vector{Int}[]
  sizehint!(coords, d)
  for g in Iterators.product([-x:x for x in max_range]...)
    r1 = collect(g)
    r2 = inverse_scale_matrix * r1
    if all(0 <= x < 1 for x in r2)
      push!(coords, r1)
    end
  end
  @assert length(coords) == d
  coord_indices = Dict(r => i for (i, r) in enumerate(coords))
  function torus_wrap(r ::AbstractVector{<:Integer})
    R = Int.(floor.(inverse_scale_matrix * r))
    r2 = r - scale_matrix * R
    return coord_indices[r2], R
  end
  @assert all(torus_wrap(r) == (i, zeros(n)) for (i, r) in enumerate(coords))
  return (coords, torus_wrap)
end
