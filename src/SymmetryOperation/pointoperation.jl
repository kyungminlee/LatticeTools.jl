export PointOperation

export apply_operation
export canonize
export iscanonical
export combinable

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
scalartype(arg::PointOperation{S}) where S = S

import Base.==
(==)(lhs::PointOperation, rhs::PointOperation) = lhs.matrix == rhs.matrix

import Base.isequal
isequal(lhs::PointOperation, rhs::PointOperation) = isequal(lhs.matrix, rhs.matrix)

import Base.*
(*)(lhs::PointOperation, rhs::PointOperation) = PointOperation(lhs.matrix * rhs.matrix)

import Base.^
function (^)(lhs::PointOperation, rhs::Integer)
    if rhs >= 0
        PointOperation(lhs.matrix^rhs)
    else
        lhs_inv_matrix = ExactLinearAlgebra.inverse(lhs.matrix)
        PointOperation(lhs_inv_matrix^(-rhs))
    end
end


combinable(lhs::PointOperation, rhs::PointOperation) = true


import Base.hash
hash(arg::PointOperation) = hash(arg.matrix)

import Base.inv
function inv(lhs::PointOperation{S}) where {S<:Real}
    PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))
end


function apply_operation(symop::PointOperation,
                        coord::AbstractVector{<:Real})
    return symop.matrix * coord
end

function (symop::PointOperation)(coord::AbstractVector{<:Real})
    return symop.matrix * coord
end


canonize(arg::PointOperation) = isone(arg.matrix) ? IdentityOperation() : arg
iscanonical(arg::PointOperation) = !isone(arg.matrix)
