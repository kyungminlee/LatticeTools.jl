export SpaceOperation
export apply_operation
export dimension
export isidentity, istranslation, ispoint

import LinearAlgebra


# S(r) = M ⋅ ( r + R )
struct SpaceOperation{Tp<:Real, Tt<:Real} <:AbstractSymmetryOperation{Tt}
    matrix::Matrix{Tp}
    displacement::Vector{Tt}

    function SpaceOperation{Tp, Tt}(dim::Integer) where {Tp<:Real, Tt<:Real}
        new{Tp, Tt}(Matrix{Tp}(LinearAlgebra.I, dim, dim), zeros(Tt, dim))
    end

    function SpaceOperation(::Type{Tp}, ::Type{Tt}, dim::Integer) where {Tp<:Real, Tt<:Real}
        new{Tp, Tt}(Matrix{Tp}(LinearAlgebra.I, dim, dim), zeros(Tt, dim))
    end

    function SpaceOperation{Tp, Tt}(matrix::AbstractMatrix{Tp},
                                    displacement::AbstractVector{Tt}) where {Tp<:Real, Tt<:Real}
        dim = size(matrix, 1)
        if size(matrix, 2) != dim
            throw(DimensionMismatch("matrix is not square: dimensions are $(size(matrix))"))
        elseif length(displacement) != dim
            throw(DimensionMismatch("dimensions of matrix and displacement do not match: matrix dimensions are $(size(matrix)), and vector dimension is $(length(displacement))"))
        end
        return new{Tp, Tt}(matrix, displacement)
    end



    function SpaceOperation(matrix::AbstractMatrix{Tp},
                            displacement::AbstractVector{Tt}) where {Tp<:Real, Tt<:Real}
        return SpaceOperation{Tp, Tt}(matrix, displacement)
    end

    function SpaceOperation(point::PointOperation{Tp},
                            translation::TranslationOperation{Tt}) where {Tp<:Real, Tt<:Real}
        return SpaceOperation{Tp, Tt}(point.matrix, translation.displacement)
    end

    function SpaceOperation(point::PointOperation{Tp}) where {Tp}
        dim = dimension(point)
        return SpaceOperation{Tp, Tp}(point.matrix, zeros(Tp, dim))
    end

    function SpaceOperation(translation::TranslationOperation{Tt}) where {Tt}
        dim = dimension(translation)
        return SpaceOperation{Int, Tt}(Matrix{Int}(LinearAlgebra.I, dim, dim), translation.displacement)
    end
end


import Base.convert
function convert(::Type{SpaceOperation{Tp, Tt}}, op::IdentityOperation{<:Union{Tp, Tt}}) where {Tp, Tt}
    return SpaceOperation(Tp, Tt, dimension(op))
end

function convert(::Type{SpaceOperation{Tp, Tt}}, op::PointOperation{Tp}) where {Tp, Tt}
    return SpaceOperation(op)
end

function convert(::Type{SpaceOperation{Tp, Tt}}, op::TranslationOperation{Tt}) where {Tp, Tt}
    return SpaceOperation(op)
end

function convert(::Type{SpaceOperation{Tp, Tt}}, arg::AbstractVector{Tt}) where {Tp, Tt}
    dim = length(arg)
    return SpaceOperation(Matrix{Int}(LinearAlgebra.I, dim, dim), arg)
end

function convert(::Type{SpaceOperation{Tp, Tt}}, arg::AbstractMatrix{Tp}) where {Tp, Tt}
    dim = size(arg, 1)
    return SpaceOperation(arg, zeros(Tp, dim))
end


import Base.promote_rule
function promote_rule(::Type{SpaceOperation{Tp, Tt}}, ::Type{<:IdentityOperation{<:Union{Tp, Tt}}}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}
end
function promote_rule(::Type{SpaceOperation{Tp, Tt}}, ::Type{TranslationOperation{Tt}}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}
end
function promote_rule(::Type{SpaceOperation{Tp, Tt}}, ::Type{PointOperation{Tp}}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}
end
function promote_rule(::Type{TranslationOperation{Tt}}, ::Type{PointOperation{Tp}}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}
end
# function promote_rule(::Type{PointOperation{Tp}}, ::Type{TranslationOperation{Tt}}) where {Tp, Tt}
#     return SpaceOperation{Tp, Tt}
# end


## properties
dimension(arg::SpaceOperation) = length(arg.displacement)
isidentity(arg::SpaceOperation) = isone(arg.matrix) && iszero(arg.displacement)
istranslation(arg::SpaceOperation) = isone(arg.matrix)
ispoint(arg::SpaceOperation) = iszero(arg.displacement)

function hash(op::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    h = hash(SpaceOperation{Tp, Tt})
    h = hash(op.matrix, h)
    h = hash(op.displacement, h)
    return h
end


## operators
import Base.==
function (==)(lhs::SpaceOperation{Tp, Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return lhs.matrix == rhs.matrix && lhs.displacement == rhs.displacement
end

function (==)(lhs::SpaceOperation{Tp, Tt}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    return isone(lhs.matrix) && lhs.displacement == rhs.displacement
end

function (==)(lhs::SpaceOperation{Tp, Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    return iszero(lhs.displacement) && lhs.matrix == rhs.matrix
end

function (==)(lhs::TranslationOperation{Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return isone(rhs.matrix) && lhs.displacement == rhs.displacement
end

function (==)(lhs::PointOperation{Tp}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return iszero(rhs.displacement) && lhs.matrix == rhs.matrix
end

function (==)(sop::SpaceOperation{Tp, Tt}, iden::IdentityOperation{<:Union{Tp, Tt}}) where {Tp, Tt}
    return dimension(sop) == dimension(iden) && iszero(sop.displacement) && isone(sop.matrix)
end

function (==)(iden::IdentityOperation{<:Union{Tp, Tt}}, sop::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return dimension(sop) == dimension(iden) && iszero(sop.displacement) && isone(sop.matrix)
end

import Base.*
function (*)(lhs::PointOperation{Tp}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}(lhs.matrix, rhs.displacement)
end

function (*)(lhs::TranslationOperation{Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    rhs_matrix_inv = Matrix{Tp}(ExactLinearAlgebra.inverse(rhs.matrix))
    return SpaceOperation{Tp, Tt}(rhs.matrix, rhs_matrix_inv * lhs.displacement)
end

#   ML ⋅ ( MR ⋅ ( x + ρR ) + ρL )
# = ML MR x + ML MR ρR + ML ρL
# = ML MR (x + ρR + inv(MR) ⋅ ρL )
function (*)(lhs::SpaceOperation{Tp, Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function (*)(lhs::SpaceOperation{Tp, Tt}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    matrix = lhs.matrix
    displacement = rhs.displacement + lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function (*)(lhs::TranslationOperation{Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function (*)(lhs::SpaceOperation{Tp, Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function (*)(lhs::PointOperation{Tp}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    displacement = rhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

import Base.inv
function inv(arg::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix_inv = Matrix{Tp}(ExactLinearAlgebra.inverse(arg.matrix))
    return SpaceOperation{Tp, Tt}(matrix_inv, -arg.matrix * arg.displacement)
end

import Base.^
function (^)(op::SpaceOperation{Tp, Tt}, power::Integer) where {Tp, Tt}
    if power == 0
        return SpaceOperation(Tp, Tt, dimension(op))
    elseif power < 0
        op_inv = inv(op)
        return op_inv^(-power)
    else
        return op * op^(power-1)
    end
end


## apply
function apply_operation(op::SpaceOperation{Tp, Tt}, coord::AbstractArray{<:Union{Tp, Tt}}) where {Tp, Tt}
    return op.matrix * (coord .+ op.displacement)
end

function (op::SpaceOperation{Tp, Tt})(coord::AbstractArray{<:Union{Tp, Tt}}) where {Tp, Tt}
    return op.matrix * (coord .+ op.displacement)
end


