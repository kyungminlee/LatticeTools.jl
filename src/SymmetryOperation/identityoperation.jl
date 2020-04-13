export IdentityOperation

export inverse
export apply_symmetry
export canonize
export iscanonical

struct IdentityOperation <:AbstractSymmetryOperation end

import Base.*
(*)(lhs::AbstractSymmetryOperation, rhs::IdentityOperation) = lhs
(*)(lhs::IdentityOperation, rhs::AbstractSymmetryOperation) = rhs


inverse(arg::IdentityOperation) = arg
apply_symmetry(symop::IdentityOperation, rhs) = rhs
(symop::IdentityOperation)(coord::AbstractVector{<:Real}) = coord

canonize(arg::IdentityOperation) = arg

iscanonical(arg::IdentityOperation) = true