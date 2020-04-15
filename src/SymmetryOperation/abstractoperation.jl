export AbstractSymmetryOperation
export combinable

abstract type AbstractSymmetryOperation{S<:Real} end

combinable(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation) = false
domaintype(arg::AbstractSymmetryOperation{S}) where {S<:Real} = S