export HypercubicLattice
export get_generators

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
    @assert all(torus_wrap(r) == (zeros(n), i) for (i, r) in enumerate(coords))

    return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, torus_wrap)
  end
end
