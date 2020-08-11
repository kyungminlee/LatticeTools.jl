export AbstractSpaceSymmetryOperation
export domaintype

"""
    AbstractSpaceSymmetryOperation{S<:Real}

Abstract space symmetry operation, i.e. translation, point, and space operation.
See also [`TranslationOperation`](@ref), [`PointOperation`](@ref), [`SpaceOperation`](@ref).
"""
abstract type AbstractSpaceSymmetryOperation{S<:Real} end

"""
    domaintype(arg::AbstractSpaceSymmetryOperation{S}) where {S<:Real}

Domain type of `arg`, i.e. the type of the coordinates.
"""
domaintype(arg::AbstractSpaceSymmetryOperation{S}) where {S<:Real} = S
