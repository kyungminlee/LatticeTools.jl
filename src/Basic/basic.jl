export ExactLinearAlgebra

module ExactLinearAlgebra
    import LinearAlgebra

    ScalarType = Union{Integer, Complex{<:Integer}, <:Rational, Complex{<:Rational}}

    function get_cofactor_matrix_unsafe!(
        out::AbstractMatrix{I1},
        mat::AbstractMatrix{I2},
        row::Integer,
        col::Integer,
    ) where {I1<:ScalarType, I2<:ScalarType}
        out[1:row-1, 1:col-1] = mat[1:row-1, 1:col-1]
        out[1:row-1, col:end] = mat[1:row-1, col+1:end]
        out[row:end, 1:col-1] = mat[row+1:end, 1:col-1]
        out[row:end, col:end] = mat[row+1:end, col+1:end]
        return out
    end

    function determinant_unsafe(mat::AbstractMatrix{S}) where {S<:ScalarType}
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

    function determinant(mat::AbstractMatrix{S}) where {S<:ScalarType}
        n, m = size(mat)
        n != m && throw(ArgumentError("matrix needs to be square"))
        n <= 0 && throw(ArgumentError("matrix empty"))
        return determinant_unsafe(mat)
    end

    determinant(mat::AbstractMatrix{<:AbstractFloat}) = LinearAlgebra.det(mat)

    function inverse(mat::AbstractMatrix{S}) where {S<:ScalarType}
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

    inverse(mat::AbstractMatrix{<:AbstractFloat}) = LinearAlgebra.inv(mat)

end # module ExactLinearAlgebra


module ExprParser

export parse_expr

SYMBOL_DATABASE = Dict(
    :i => im,
    :im => im,
    :pi => pi,
    :π => π,
    :+ => (+),  :- => (-),  :* => (*),  :/ => (/),  :\ => (\),  :^ => (^), :// => (//),
    :cos => cos,  :sin => sin,  :tan => tan,
    :exp => exp,  :cis => cis,
    :cospi => cospi,  :sinpi => sinpi,
    :sqrt => sqrt,
    :log => log,
    :abs => abs,  :abs2 => abs2,  :sign => sign,
    :conj => conj, :real => real, :imag => imag,
    :angle => angle
)

myeval(expr::Number) = expr
function myeval(expr::Symbol)
    if haskey(SYMBOL_DATABASE, expr)
        return SYMBOL_DATABASE[expr]
    else
        error("unsupported symbol $expr")
    end
end

function myeval(expr::Expr)
    if expr.head == :call
        ftn, args = Iterators.peel(myeval.(expr.args))
        return ftn(args...)
    elseif expr.head == :vect
        return [(myeval.(expr.args))...]
    elseif expr.head == :tuple
        return tuple((myeval.(expr.args))...)
    elseif expr.head == :hcat
        return hcat((myeval.(expr.args))...)
    elseif expr.head == :row
        return hcat((myeval.(expr.args))...)
    elseif expr.head == :vcat
        return vcat((myeval.(expr.args))...)
    else
        error("unsupported expression $expr")
    end
end

parse_expr(expr::Number) = expr
parse_expr(expr::AbstractString) = myeval(Meta.parse(expr))
parse_expr(expr::AbstractArray) = [parse_expr(elem) for elem in expr]

end


using .ExprParser: parse_expr

#= # commented for now
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
=#


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
cleanup_number(x::Rational{<:Integer}, tol::Real) = x

function cleanup_number(x::Complex, tol::Real)
    return cleanup_number(real(x), tol) + im * cleanup_number(imag(x), tol)
end

cleanup_number(x::AbstractArray, tol::Real) = [cleanup_number(y, tol) for y in x]


function extended_gcd(a::Integer, b::Integer)
    sign_a, sign_b = sign(a), sign(b)
    s, old_s = 0, 1
    r, old_r = abs(b), abs(a)
    while r != 0
        quotient = old_r ÷ r
        (old_r, r) = (r, old_r - quotient * r)
        (old_s, s) = (s, old_s - quotient * s)
    end
    if b != 0
        bezout_t = (old_r - old_s * abs(a)) ÷ abs(b)
    else
        bezout_t = 0
    end
    return (old_r, [sign_a * old_s, sign_b * bezout_t])
end
