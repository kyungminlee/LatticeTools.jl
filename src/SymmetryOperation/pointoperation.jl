export PointOperation

export apply_operation
export domaintype
export isidentity, istranslation, ispoint

"""
    PointOperation{S<:Real}

Represents point symmetry operation

# Fields
* `matrix`: rotation/mirror/inversion matrix
"""
struct PointOperation{S<:Real} <:AbstractSpaceSymmetryOperation{S}
    matrix::Matrix{S}
    function PointOperation(matrix::AbstractMatrix{S}) where {S<:Real}
        return PointOperation{S}(matrix)
    end
    function PointOperation{S}(matrix::AbstractMatrix{S2}) where {S<:Real, S2<:Real}
        dim = size(matrix, 1)
        size(matrix, 2) != dim && throw(DimensionMismatch("matrix not square"))
        return new{S}(matrix)
    end
end


function Base.convert(::Type{PointOperation{S}}, obj::IdentityOperation{S}) where S
    dim = dimension(obj)
    return PointOperation{S}(Matrix(LinearAlgebra.I, dim, dim))
end

function Base.convert(::Type{PointOperation{S}}, matrix::AbstractMatrix{S}) where S
    return PointOperation{S}(matrix)
end


function Base.promote_rule(::Type{PointOperation{S}}, ::Type{IdentityOperation{S}}) where S
    return PointOperation{S}
end


## properties
dimension(arg::PointOperation) = size(arg.matrix, 1)
isidentity(arg::PointOperation) = isone(arg.matrix)
istranslation(arg::PointOperation) = isone(arg.matrix)
ispoint(arg::PointOperation) = true


Base.hash(arg::PointOperation{S}) where S = hash(arg.matrix, hash(PointOperation{S}))


## operators
function Base.:(==)(lhs::PointOperation{S}, rhs::PointOperation{S}) where S
    return lhs.matrix == rhs.matrix
end
function Base.:(==)(pop::PointOperation{S}, iden::IdentityOperation{S}) where S
    return isidentity(pop)
end
function Base.:(==)(iden::IdentityOperation{S}, pop::PointOperation{S}) where S
    return isidentity(pop)
end
function Base.:(==)(pop::PointOperation{S}, top::TranslationOperation{S}) where S
    return isidentity(pop) && isidentity(top)
end
function Base.:(==)(top::TranslationOperation{S}, pop::PointOperation{S}) where S
    return isidentity(pop) && isidentity(top)
end

function Base.:(*)(lhs::PointOperation, rhs::PointOperation)
    return PointOperation(lhs.matrix * rhs.matrix)
end

function Base.:(^)(lhs::PointOperation{S}, rhs::Integer) where S
    if rhs >= 0
        return PointOperation(lhs.matrix^rhs)
    else
        lhs_inv_matrix = Matrix{S}(ExactLinearAlgebra.inverse(lhs.matrix))
        return PointOperation{S}(lhs_inv_matrix^(-rhs))
    end
end

function Base.inv(lhs::PointOperation{S}) where {S<:Real}
    return PointOperation{S}(ExactLinearAlgebra.inverse(lhs.matrix))
end


## apply
"""
    apply_operation(symop::PointOperation{S}, coordinate::AbstractArray{S}) where {S}

Apply point operation to the coordinates.
"""
function apply_operation(symop::PointOperation{S}, coord::AbstractVector{S}) where S
    return symop.matrix * coord
end

function (symop::PointOperation{S})(coord::AbstractVector{S}) where S
    return symop.matrix * coord
end
