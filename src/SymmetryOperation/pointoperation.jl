export PointOperation

export apply_operation
export canonize
export iscanonical
export combinable
export domaintype
export isidentity

struct PointOperation{S<:Real} <:AbstractSymmetryOperation{S}
    matrix::Matrix{S}

    function PointOperation(matrix::AbstractMatrix{S}) where {S<:Real}
        return PointOperation{S}(matrix)
    end

    function PointOperation{S}(matrix::AbstractMatrix{S2}) where {S<:Real, S2<:Real}
        dim = size(matrix, 1)
        size(matrix, 2) != dim && throw(DimensionMismatch("matrix not square"))
        new{S}(matrix)
    end
end


## properties
dimension(arg::PointOperation) = size(arg.matrix, 1)
# domaintype(arg::PointOperation{S}) where S = S
combinable(lhs::PointOperation{S}, rhs::PointOperation{S}) where S = dimension(lhs) == dimension(rhs)
isidentity(arg::PointOperation) = isone(arg.matrix)

import Base.hash
hash(arg::PointOperation) = hash(arg.matrix)


## operators
import Base.==
(==)(lhs::PointOperation, rhs::PointOperation) = lhs.matrix == rhs.matrix
(==)(pop::PointOperation, iden::IdentityOperation) = isone(pop.matrix)
(==)(iden::IdentityOperation, pop::PointOperation) = isone(pop.matrix)
(==)(pop::PointOperation, top::TranslationOperation) = isone(pop.matrix) && iszero(top.displacement)
(==)(top::TranslationOperation, pop::PointOperation) = isone(pop.matrix) && iszero(top.displacement)

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

import Base.inv
function inv(lhs::PointOperation{S}) where {S<:Real}
    PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))
end

## apply
apply_operation(symop::PointOperation, coord::AbstractVector{<:Real}) = symop.matrix * coord
(symop::PointOperation)(coord::AbstractVector{<:Real}) = symop.matrix * coord

# ## canonical
# canonize(arg::PointOperation) = isone(arg.matrix) ? IdentityOperation() : arg
# iscanonical(arg::PointOperation) = !isone(arg.matrix)

