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
(==)(lhs::PointOperation, rhs::PointOperation) = lhs.matrix == rhs.matrix

import Base.*
(*)(lhs::PointOperation, rhs::PointOperation) = PointOperation(lhs.matrix * rhs.matrix)
(*)(lhs::PointOperation, rhs::Real) = PointOperation(lhs.matrix .* rhs)
(*)(lhs::Real, rhs::PointOperation) = PointOperation(lhs .* rhs.matrix)


inverse(lhs::PointOperation{S}) where {S<:Real} = PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))


function apply_symmetry(symop::PointOperation,
                        coord::AbstractVector{<:Real})
    return symop.matrix * coord
end

function (symop::PointOperation)(coord::AbstractVector{<:Real})
    return symop.matrix * coord
end


canonize(arg::PointOperation) = isone(arg.matrix) ? IdentityOperation() : arg
iscanonical(arg::PointOperation) = !isone(arg.matrix)