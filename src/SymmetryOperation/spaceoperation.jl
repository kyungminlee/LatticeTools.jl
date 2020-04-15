export SpaceOperation
export apply_operation
export canonize
export iscanonical

using LinearAlgebra


# S(r) = M ⋅ ( r + R )
struct SpaceOperation{S<:Real} <:AbstractSymmetryOperation{S}
    matrix::Matrix{S}
    displacement::Vector{S}

    function SpaceOperation{S}(dim::Integer) where {S<:Real}
        new{S}(Matrix(I, dim, dim), zeros(Bool, dim))
    end

    function SpaceOperation(::Type{S}, dim::Integer) where {S<:Real}
        new{S}(Matrix(I, dim, dim), zeros(Bool, dim))
    end

    function SpaceOperation(matrix::AbstractMatrix{S},
                            displacement::AbstractVector{S}) where {S<:Real}
        dim = size(matrix, 1)
        if size(matrix, 2) != dim
            throw(DimensionMismatch("matrix is not square: dimensions are $(size(matrix))"))
        elseif length(displacement) != dim
            throw(DimensionMismatch("dimensions of matrix and displacement do not match: matrix dimensions are $(size(matrix)), and vector dimension is $(length(displacement))"))
        end
        return new{S}(matrix, displacement)
    end

    function SpaceOperation(point::PointOperation{S},
                            translation::TranslationOperation{S}) where {S<:Real}
        return new{S}(point.matrix, translation.displacement)
    end

    function SpaceOperation(point::PointOperation{S}) where S
        dim = dimension(point)
        return new{S}(point.matrix, zeros(Int, dim))
    end

    function SpaceOperation(translation::TranslationOperation{S}) where S
        dim = dimension(translation)
        return new{S}(Matrix(I, dim, dim), translation.displacement)
    end

    # function SpaceOperation(arg::ProductOperation{S, F}) where {S, F}
    #     dim = dimension(arg)
    #     if dim == 0
    #         throw(ArgumentError("Cannot create space operation out of an identity (dimension undetermined)"))
    #     end
        
    #     canon = canonize(arg)
    #     if isa(canon, IdentityOperation)
    #         return SpaceOperation(Matrix(I, dim, dim), zeros(S, dim))
    #     elseif isa(canon, TranslationOperation)
    #         return SpaceOperation(Matrix(I, dim, dim), canon.displacement)
    #     elseif isa(canon, PointOperation)
    #         return SpaceOperation(canon.matrix, zeros(S, dim))
    #     else
    #         @assert length(canon.factors) == 2
    #         return SpaceOperation(canon.factors[1], canon.factors[2])
    #     end
    # end
end


import Base.convert

function convert(::Type{SpaceOperation{S}}, op::IdentityOperation{S}) where S
    return SpaceOperation{Int}(dimension(op))
end

function convert(::Type{SpaceOperation{S}}, op::PointOperation{S}) where S
    return SpaceOperation(op)
end

function convert(::Type{SpaceOperation{S}}, op::TranslationOperation{S}) where S
    return SpaceOperation(op)
end



import Base.==
function (==)(lhs::SpaceOperation{S}, rhs::SpaceOperation{S}) where S
    return lhs.matrix == rhs.matrix && lhs.displacement == rhs.displacement
end

function (==)(lhs::SpaceOperation{S}, rhs::TranslationOperation{S}) where S
    return isone(lhs.matrix) && lhs.displacement == rhs.displacement
end

function (==)(lhs::SpaceOperation{S}, rhs::PointOperation{S}) where S
    return iszero(lhs.displacement) && lhs.matrix == rhs.matrix
end

function (==)(lhs::TranslationOperation{S}, rhs::SpaceOperation{S}) where S
    return isone(rhs.matrix) && lhs.displacement == rhs.displacement
end

function (==)(lhs::PointOperation{S}, rhs::SpaceOperation{S}) where S
    return iszero(rhs.displacement) && lhs.matrix == rhs.matrix
end

function (==)(sop::SpaceOperation{S}, iden::IdentityOperation) where S
    return iszero(sop.displacement) && isone(sop.matrix)
end

function (==)(iden::IdentityOperation, sop::SpaceOperation{S}) where S
    return iszero(sop.displacement) && isone(sop.matrix)
end


import Base.*
function (*)(lhs::PointOperation{S}, rhs::TranslationOperation{S}) where S
    return SpaceOperation(lhs.matrix, rhs.displacement)
end

function (*)(lhs::TranslationOperation{S}, rhs::PointOperation{S}) where S
    rhs_matrix_inv = Matrix{S}(ExactLinearAlgebra.inverse(rhs.matrix))
    return SpaceOperation(rhs.matrix, rhs_matrix_inv * lhs.displacement)
end

#   ML ⋅ ( MR ⋅ ( x + ρR ) + ρL )
# = ML MR x + ML MR ρR + ML ρL
# = ML MR (x + ρR + inv(MR) ⋅ ρL )
function (*)(lhs::SpaceOperation{S}, rhs::SpaceOperation{S}) where {S}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{S}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation(matrix, displacement)
end

function (*)(lhs::SpaceOperation{S}, rhs::TranslationOperation{S}) where {S}
    matrix = lhs.matrix
    displacement = rhs.displacement + lhs.displacement
    return SpaceOperation(matrix, displacement)
end

function (*)(lhs::TranslationOperation{S}, rhs::SpaceOperation{S}) where {S}
    matrix = rhs.matrix
    rhs_matrix_inv = Matrix{S}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation(matrix, displacement)
end

function (*)(lhs::SpaceOperation{S}, rhs::PointOperation{S}) where {S}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{S}(ExactLinearAlgebra.inverse(rhs.matrix))
    displacement = rhs_matrix_inv * lhs.displacement
    return SpaceOperation(matrix, displacement)
end

function (*)(lhs::PointOperation{S}, rhs::SpaceOperation{S}) where {S}
    matrix = lhs.matrix * rhs.matrix
    displacement = rhs.displacement
    return SpaceOperation(matrix, displacement)
end

import Base.inv
function inv(arg::SpaceOperation{S}) where S
    matrix_inv = Matrix{S}(ExactLinearAlgebra.inverse(arg.matrix))
    return SpaceOperation(matrix_inv, -arg.matrix * arg.displacement)
end

