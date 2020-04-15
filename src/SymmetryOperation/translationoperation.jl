export TranslationOperation

export apply_operation
export canonize
export iscanonical
export combinable
export domaintype
export isidentity

struct TranslationOperation{S<:Real} <:AbstractSymmetryOperation{S}
    displacement::Vector{S}
    function TranslationOperation{S}(displacement::AbstractVector{<:Real}) where {S}
        new{S}(displacement)
    end
    function TranslationOperation(displacement::AbstractVector{S}) where {S<:Real}
        new{S}(displacement)
    end
end

dimension(arg::TranslationOperation) = length(arg.displacement)
# domaintype(arg::TranslationOperation{S}) where S = S

import Base.==
function (==)(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S
    lhs.displacement == rhs.displacement
end
function (==)(top::TranslationOperation{S}, iden::IdentityOperation{S}) where S
    iszero(top.displacement)
end
function (==)(iden::IdentityOperation{S}, top::TranslationOperation{S}) where S
    iszero(top.displacement)
end

import Base.isequal
function isequal(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S
    isequal(lhs.displacement, rhs.displacement)
end

import Base.*
function (*)(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S
    TranslationOperation{S}(lhs.displacement .+ rhs.displacement)
end

import Base.^
function (^)(lhs::TranslationOperation{S}, rhs::Real) where S
    TranslationOperation{S}(lhs.displacement .* rhs)
end

import Base.inv
function inv(arg::TranslationOperation{S}) where S
    TranslationOperation{S}(-arg.displacement)
end

combinable(lhs::TranslationOperation{S}, rhs::TranslationOperation{S}) where S = true


import Base.hash
hash(arg::TranslationOperation) = hash(arg.displacement)


function apply_operation(symop::TranslationOperation{S},
                         coord::AbstractArray{S}) where {S<:Real}
    return coord .+ symop.displacement
end

function (symop::TranslationOperation{S})(coord::AbstractVector{S}) where S
    return coord .+ symop.displacement
end


# function canonize(arg::TranslationOperation{S}) where S
#     if iszero(arg.displacement) 
#         IdentityOperation{S}(dimension(arg))
#     else
#         arg
#     end
# end

# iscanonical(arg::TranslationOperation) = !iszero(arg.displacement)

isidentity(arg::TranslationOperation) = iszero(arg.displacement)