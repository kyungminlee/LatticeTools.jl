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

    @doc """
        IdentityOperation{S}(dim::Integer) where {S<:Real}

    Construct an identity operation of dimension `dim`, on coordinates of type `S`.

    * 2020-08-10: Haven't figured out how to make this appear in the documentation.
    """
    function IdentityOperation{S}(dim::Integer) where {S<:Real}
        return new{S}(dim)
    end

    @doc """
        IdentityOperation(S, dim::Integer)

    Construct an identity operation of dimension `dim`, on coordinates of type `S`.
    """
    function IdentityOperation(::Type{S}, dim::Integer) where {S<:Real}
        return new{S}(dim)
    end
end


## properties

"""
    isidentity(arg::IdentityOperation)

Check whether the argument is an identity. Always `true`.
"""
isidentity(arg::IdentityOperation) = true


"""
    istranslation(arg::IdentityOperation)

Check whether the argument is a translation operation. Always `true`.
"""
istranslation(arg::IdentityOperation) = true


"""
    ispoint(arg::IdentityOperation)

Check whether the argument is a point operation. Always `true`.
"""
ispoint(arg::IdentityOperation) = true

"""
    dimension(arg::IdentityOperation)

Return the spatial dimension of the identity operation.
"""
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
