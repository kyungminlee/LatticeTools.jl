export TranslationOperation

export apply_operation
export domaintype
export isidentity, istranslation, ispoint

struct TranslationOperation{S<:Real} <:AbstractSymmetryOperation{S}
    displacement::Vector{S}
    function TranslationOperation{S}(displacement::AbstractVector{<:Real}) where {S}
        new{S}(displacement)
    end
    function TranslationOperation(displacement::AbstractVector{S}) where {S<:Real}
        new{S}(displacement)
    end
end


function Base.convert(::Type{TranslationOperation{S}}, obj::IdentityOperation{S}) where S
    dim = dimension(obj)
    return TranslationOperation{S}(zeros(S, dim))
end

function Base.convert(::Type{TranslationOperation{S}}, displacement::AbstractVector{S}) where S
    return TranslationOperation{S}(displacement)
end


function Base.promote_rule(::Type{TranslationOperation{S}}, ::Type{IdentityOperation{S}}) where S
    return TranslationOperation{S}
end


## properties
dimension(arg::TranslationOperation) = length(arg.displacement)
isidentity(arg::TranslationOperation) = iszero(arg.displacement)
istranslation(arg::TranslationOperation) = true
ispoint(arg::TranslationOperation) = iszero(arg.displacement)


Base.hash(arg::TranslationOperation{S}) where S = hash(arg.displacement, hash(TranslationOperation{S}))


## operators
function Base.:(==)(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S
    lhs.displacement == rhs.displacement
end
function Base.:(==)(top::TranslationOperation{S}, iden::IdentityOperation{S}) where S
    iszero(top.displacement)
end
function Base.:(==)(iden::IdentityOperation{S}, top::TranslationOperation{S}) where S
    iszero(top.displacement)
end


function Base.:(*)(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S
    TranslationOperation{S}(lhs.displacement .+ rhs.displacement)
end


function Base.:(^)(lhs::TranslationOperation{S}, rhs::Real) where S
    TranslationOperation{S}(lhs.displacement .* rhs)
end


function Base.inv(arg::TranslationOperation{S}) where S
    TranslationOperation{S}(-arg.displacement)
end


## apply
"""
    apply_operation(symop::TranslationOperation{S}, coordinate::AbstractArray{S}) where S
"""
function apply_operation(symop::TranslationOperation{S},
                         coord::AbstractArray{S}) where {S<:Real}
    return coord .+ symop.displacement
end

function (symop::TranslationOperation{S})(coord::AbstractVector{S}) where S
    return coord .+ symop.displacement
end
