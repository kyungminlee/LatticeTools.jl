export PointOperation

export apply_operation
export domaintype
export isidentity, istranslation, ispoint

import LinearAlgebraX

"""
    PointOperation{S<:Real}

Represents point symmetry operation

# Fields
* `matrix`: rotation/mirror/inversion matrix
"""
struct PointOperation{S<:Real} <:AbstractSpaceSymmetryOperation{S}
    matrix::Matrix{S}

    @doc """
        PointOperation(matrix::AbstractMatrix{S}) where {S<:Real}

    Construct a point operation with `matrix`.
    """
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

"""
    isidentity(arg::PointOperation)

Check whether `arg` is identity.
"""
isidentity(arg::PointOperation) = isone(arg.matrix)

"""
    istranslation(arg::PointOperation)

Check whether `arg` is a translation operation, i.e. identity.
"""
istranslation(arg::PointOperation) = isone(arg.matrix)

"""
    ispoint(arg::PointOperation)

Check whether `arg` is a point operation. Always `true`.
"""
ispoint(arg::PointOperation) = true

"""
    dimension(arg::PointOperation)

Return spatial dimension of `arg`.
"""
dimension(arg::PointOperation) = size(arg.matrix, 1)

Base.hash(arg::PointOperation{S}, h::UInt) where S = hash(PointOperation{S}, hash(arg.matrix, h))


## operators
function Base.:(==)(lhs::PointOperation, rhs::PointOperation)
    return lhs.matrix == rhs.matrix
end
function Base.:(==)(pop::PointOperation, ::IdentityOperation)
    return isidentity(pop)
end
function Base.:(==)(::IdentityOperation, pop::PointOperation)
    return isidentity(pop)
end
function Base.:(==)(pop::PointOperation, top::TranslationOperation)
    return isidentity(pop) && isidentity(top)
end
function Base.:(==)(top::TranslationOperation, pop::PointOperation)
    return isidentity(pop) && isidentity(top)
end

function Base.:(*)(lhs::PointOperation, rhs::PointOperation)
    return PointOperation(lhs.matrix * rhs.matrix)
end

function Base.:(^)(lhs::PointOperation{S}, rhs::Integer) where S
    if rhs >= 0
        return PointOperation(lhs.matrix^rhs)
    else
        lhs_inv_matrix = Matrix{S}(LinearAlgebraX.invx(lhs.matrix))
        return PointOperation{S}(lhs_inv_matrix^(-rhs))
    end
end

function Base.inv(lhs::PointOperation{S}) where {S<:Real}
    return PointOperation{S}(LinearAlgebraX.invx(lhs.matrix))
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
