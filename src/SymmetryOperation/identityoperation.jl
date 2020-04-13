export IdentityOperation

export apply_operation
export canonize
export iscanonical
export combinable

struct IdentityOperation <:AbstractSymmetryOperation end

import Base.*
(*)(lhs::AbstractSymmetryOperation, rhs::IdentityOperation) = lhs
(*)(lhs::IdentityOperation, rhs::AbstractSymmetryOperation) = rhs

combinable(lhs::AbstractSymmetryOperation, rhs::IdentityOperation) = true
combinable(lhs::IdentityOperation, rhs::AbstractSymmetryOperation) = true

import Base.inv
inv(arg::IdentityOperation) = arg
apply_operation(symop::IdentityOperation, rhs) = rhs
(symop::IdentityOperation)(coord::AbstractVector{<:Real}) = coord

canonize(arg::IdentityOperation) = arg

iscanonical(arg::IdentityOperation) = true