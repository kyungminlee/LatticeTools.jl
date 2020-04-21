export AbstractSymmetryOperation
export domaintype

abstract type AbstractSymmetryOperation{S<:Real} end

domaintype(arg::AbstractSymmetryOperation{S}) where {S<:Real} = S