export IdentityOperation

export apply_operation
export domaintype
export isidentity, istranslation, ispoint

"""
    IdentityOperation{S<:Real} <: AbstractSpaceSymmetryOperation{S}

Represents identity (space symmetry) operation

# Fields
* `dimension::Int`: dimension of the space on which the identity operation acts
"""
struct IdentityOperation{S<:Real} <: AbstractSpaceSymmetryOperation{S}
    dimension::Int
    function IdentityOperation{S}(dim::Integer) where {S<:Real}
        return new{S}(dim)
    end
    function IdentityOperation(::Type{S}, dim::Integer) where {S<:Real}
        return new{S}(dim)
    end
end


## properties
isidentity(arg::IdentityOperation) = true
istranslation(arg::IdentityOperation) = true
ispoint(arg::IdentityOperation) = true

dimension(arg::IdentityOperation) = arg.dimension


function Base.hash(arg::IdentityOperation{S}) where S
    return hash(arg.dimension, hash(IdentityOperation{S}))
end


## operators
function Base.:(==)(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where S
    return lhs.dimension == rhs.dimension
end

function Base.:(*)(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where S
    if dimension(lhs) != dimension(rhs)
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return lhs
end

function Base.:(*)(lhs::AbstractSpaceSymmetryOperation{S}, rhs::IdentityOperation{S}) where S
    if dimension(lhs) != dimension(rhs)
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return lhs
end

function Base.:(*)(lhs::IdentityOperation{S}, rhs::AbstractSpaceSymmetryOperation{S}) where S
    if dimension(lhs) != dimension(rhs)
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return rhs
end


Base.:(^)(lhs::IdentityOperation, rhs::Integer) = lhs


Base.inv(arg::IdentityOperation) = arg


## apply
"""
    apply_operation(identity{S}, coordinate::AbstractArray{S}) where {S<:Real}

Do nothing.
"""
function apply_operation(symop::IdentityOperation{S}, arg::AbstractArray{S}) where {S<:Real}
    if dimension(symop) != size(arg, 1)
        throw(DimensionMismatch("dimension mismatch"))
    end
    return arg
end

function (symop::IdentityOperation{S})(arg::AbstractArray{S}) where {S<:Real}
    if dimension(symop) != size(arg, 1)
        throw(DimensionMismatch("dimension mismatch"))
    end
    return arg
end
