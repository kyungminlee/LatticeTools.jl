export IdentityOperation

export apply_operation
export canonize
export iscanonical
export combinable
export domaintype
export isidentity

struct IdentityOperation{S<:Real} <: AbstractSymmetryOperation{S}
    dimension::Int
    function IdentityOperation{S}(dim::Integer) where {S<:Real} 
        new{S}(dim)
    end
    function IdentityOperation(::Type{S}, dim::Integer) where {S<:Real} 
        new{S}(dim)
    end
end

import Base.==
function (==)(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where S 
    lhs.dimension == rhs.dimension
end

import Base.*
function (*)(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where S
    if dimension(lhs) != dimension(rhs) 
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return lhs
end

function (*)(lhs::AbstractSymmetryOperation{S}, rhs::IdentityOperation{S}) where S
    if dimension(lhs) != dimension(rhs) 
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return lhs
end

function (*)(lhs::IdentityOperation{S}, rhs::AbstractSymmetryOperation{S}) where S
    if dimension(lhs) != dimension(rhs) 
        throw(DimensionMismatch("dimensions mismatch"))
    end
    return rhs
end

import Base.^
(^)(lhs::IdentityOperation, rhs::Integer) = lhs

import Base.inv
inv(arg::IdentityOperation) = arg

# combinable(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where {S<:Real} = dimension(lhs) == dimension(rhs)
# combinable(lhs::AbstractSymmetryOperation{S}, rhs::IdentityOperation{S}) where {S<:Real} = dimension(lhs) == dimension(rhs)
# combinable(lhs::IdentityOperation{S}, rhs::AbstractSymmetryOperation{S}) where {S<:Real} = dimension(lhs) == dimension(rhs)

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

isidentity(arg::IdentityOperation) = true

# canonize(arg::IdentityOperation) = arg
# iscanonical(arg::IdentityOperation) = true

dimension(arg::IdentityOperation) = arg.dimension
domaintype(arg::IdentityOperation{S}) where {S<:Real} = S
