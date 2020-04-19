export PointOperation

export apply_operation
export domaintype
export isidentity, istranslation, ispoint

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


import Base.convert
function convert(::Type{PointOperation{S}}, obj::IdentityOperation{S}) where S
    dim = dimension(obj)
    return PointOperation{S}(Matrix(I, dim, dim))
end

function convert(::Type{PointOperation{S}}, matrix::AbstractMatrix{S}) where S
    return PointOperation{S}(matrix)
end


import Base.promote_rule
function promote_rule(::Type{PointOperation{S}}, ::Type{IdentityOperation{S}}) where S
    return PointOperation{S}
end


## properties
dimension(arg::PointOperation) = size(arg.matrix, 1)
isidentity(arg::PointOperation) = isone(arg.matrix)
istranslation(op::PointOperation) = isone(arg.matrix)
ispoint(op::PointOperation) = true


import Base.hash
hash(arg::PointOperation{S}) where S = hash(arg.matrix, hash(PointOperation{S}))


## operators
import Base.==
function (==)(lhs::PointOperation{S}, rhs::PointOperation{S}) where S
    lhs.matrix == rhs.matrix
end
function (==)(pop::PointOperation{S}, iden::IdentityOperation{S}) where S
    isidentity(pop)
end
function (==)(iden::IdentityOperation{S}, pop::PointOperation{S}) where S
    isidentity(pop)
end
function (==)(pop::PointOperation{S}, top::TranslationOperation{S}) where S
    isidentity(pop) && isidentity(top)
end
function (==)(top::TranslationOperation{S}, pop::PointOperation{S}) where S
    isidentity(pop) && isidentity(top)
end

import Base.isequal
isequal(lhs::PointOperation, rhs::PointOperation) = isequal(lhs.matrix, rhs.matrix)

import Base.*
(*)(lhs::PointOperation, rhs::PointOperation) = PointOperation(lhs.matrix * rhs.matrix)

import Base.^
function (^)(lhs::PointOperation{S}, rhs::Integer) where S
    if rhs >= 0
        PointOperation(lhs.matrix^rhs)
    else
        lhs_inv_matrix = Matrix{S}(ExactLinearAlgebra.inverse(lhs.matrix))
        PointOperation{S}(lhs_inv_matrix^(-rhs))
    end
end

import Base.inv
function inv(lhs::PointOperation{S}) where {S<:Real}
    PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))
end


## apply
function apply_operation(symop::PointOperation{S}, coord::AbstractVector{S}) where S
    symop.matrix * coord
end

function (symop::PointOperation{S})(coord::AbstractVector{S}) where S
    symop.matrix * coord
end

