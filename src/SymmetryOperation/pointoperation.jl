export PointOperation

export inverse
export apply_symmetry
export canonize
export iscanonical

struct PointOperation{S<:Real} <:AbstractSymmetryOperation
    matrix::Matrix{S}

    function PointOperation(matrix::AbstractMatrix{S}) where {S<:Real}
        return PointOperation{S}(matrix)
    end

    function PointOperation{S}(matrix::AbstractMatrix{S2}) where {S<:Real, S2<:Real}
        dim = size(matrix, 1)
        size(matrix, 2) != dim && throw(ArgumentError("matrix not square"))
        new{S}(matrix)
    end
end

dimension(arg::PointOperation) = size(arg.matrix, 1)

import Base.==
function (==)(lhs::PointOperation, rhs::PointOperation)
    return lhs.matrix == rhs.matrix
end

function (*)(lhs::PointOperation, rhs::PointOperation)
    return PointOperation(lhs.matrix * rhs.matrix)
end

function (*)(lhs::PointOperation, rhs::Real)
    return PointOperation(lhs.matrix .* rhs)
end

function (*)(lhs::Real, rhs::PointOperation)
    return PointOperation(lhs .* rhs.matrix)
end

inverse(lhs::PointOperation{S}) where {S<:Real} = PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))


function apply_symmetry(symop::PointOperation{S}, coord::AbstractVector{S}) where {S<:Real}
    return symop.matrix * coord
end

canonize(arg::PointOperation) = isone(arg.matrix) ? IdentityOperation() : arg
iscanonical(arg::PointOperation) = !isone(arg.matrix)