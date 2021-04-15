#export AbstractSymmetryOperation
using GroupTools

export AbstractSpaceSymmetryOperation
export domaintype

# abstract type AbstractSymmetryOperation end

"""
    AbstractSpaceSymmetryOperation{S<:Real}

Abstract space symmetry operation, i.e. translation, point, and space operation.
See also [`TranslationOperation`](@ref), [`PointOperation`](@ref), [`SpaceOperation`](@ref).
"""
abstract type AbstractSpaceSymmetryOperation{S<:Real} <: AbstractSymmetryOperation end

"""
    domaintype(arg::AbstractSpaceSymmetryOperation{S}) where {S<:Real}

Domain type of `arg`, i.e. the type of the coordinates.
"""
domaintype(::AbstractSpaceSymmetryOperation{S}) where {S<:Real} = S
domaintype(::Type{<:AbstractSpaceSymmetryOperation{S}}) where {S<:Real} = S
