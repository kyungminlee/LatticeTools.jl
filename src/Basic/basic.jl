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

parse_expr(expr::Number) = expr
function parse_expr(expr::AbstractString)
  expr2 = Meta.parse(expr)
  eval(quote
    i = im
    $expr2
  end)
end
parse_expr(expr::AbstractArray) = [parse_expr(elem) for elem in expr]

function parse_table(s::AbstractString)
  s = strip(s)
  rows = Vector{Int}[]
  for line in split(s, "\n")
    line = strip(line)
    row = [parse(Int, x) for x in split(line)]
    !isempty(row) && push!(rows, row)
  end
  return transpose(hcat(rows...))
end


function cleanup_number(x::AbstractFloat, tol::Real)
  units = [16, 16/sqrt(3), 16*sqrt(3)]
  for unit in units
    if isapprox(x*unit, round(x*unit); atol=tol)
      return round(x*unit)/unit
    end
  end
  return x
end

cleanup_number(x::Integer, tol::Real) = x
cleanup_number(x::Complex, tol::Real) = cleanup_number(real(x), tol) + im * cleanup_number(imag(x), tol)
cleanup_number(x::AbstractArray, tol::Real) = [cleanup_number(y, tol) for y in x]
