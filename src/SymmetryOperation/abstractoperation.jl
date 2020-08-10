export AbstractSpaceSymmetryOperation
export domaintype

abstract type AbstractSpaceSymmetryOperation{S<:Real} end

domaintype(arg::AbstractSpaceSymmetryOperation{S}) where {S<:Real} = S
