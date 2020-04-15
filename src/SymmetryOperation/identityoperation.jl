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
(*)(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where S = lhs
(*)(lhs::AbstractSymmetryOperation{S}, rhs::IdentityOperation{S}) where S = lhs
(*)(lhs::IdentityOperation{S}, rhs::AbstractSymmetryOperation{S}) where S = rhs

import Base.^
(^)(lhs::IdentityOperation, rhs::Integer) = lhs

import Base.inv
inv(arg::IdentityOperation) = arg

combinable(lhs::IdentityOperation{S}, rhs::IdentityOperation{S}) where {S<:Real} = true
combinable(lhs::AbstractSymmetryOperation{S}, rhs::IdentityOperation{S}) where {S<:Real} = true
combinable(lhs::IdentityOperation{S}, rhs::AbstractSymmetryOperation{S}) where {S<:Real} = true

apply_operation(symop::IdentityOperation{S}, rhs::AbstractArray{S}) where {S<:Real} = rhs
(symop::IdentityOperation{S})(arg::AbstractArray{S}) where {S<:Real} = arg

isidentity(arg::IdentityOperation) = true

# canonize(arg::IdentityOperation) = arg
# iscanonical(arg::IdentityOperation) = true

dimension(arg::IdentityOperation) = arg.dimension
domaintype(arg::IdentityOperation{S}) where {S<:Real} = S

