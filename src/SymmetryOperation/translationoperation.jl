export TranslationOperation

export inverse
export apply_symmetry
export canonize
export iscanonical
export combinable

struct TranslationOperation{S<:Real} <:AbstractSymmetryOperation
    displacement::Vector{S}
    function TranslationOperation{S}(displacement::AbstractVector{<:Real}) where {S}
        new{S}(displacement)
    end
    function TranslationOperation(displacement::AbstractVector{S}) where {S<:Real}
        new{S}(displacement)
    end
end

dimension(arg::TranslationOperation) = length(arg.displacement)

import Base.==
(==)(lhs::TranslationOperation, rhs::TranslationOperation) = lhs.displacement == rhs.displacement

import Base.isequal
isequal(lhs::TranslationOperation, rhs::TranslationOperation) = isequal(lhs.displacement, rhs.displacement)

import Base.*
(*)(lhs::TranslationOperation, rhs::TranslationOperation) = TranslationOperation(lhs.displacement .+ rhs.displacement)

import Base.^
(^)(lhs::TranslationOperation, rhs::Real) = TranslationOperation(lhs.displacement .* rhs)


combinable(lhs::TranslationOperation, rhs::TranslationOperation) = true


import Base.hash
hash(arg::TranslationOperation) = hash(arg.displacement)


inverse(arg::TranslationOperation) = TranslationOperation(-arg.displacement)


function apply_symmetry(symop::TranslationOperation, coord::AbstractVector{<:Real})
    return coord + symop.displacement
end

function (symop::TranslationOperation)(coord::AbstractVector)
    return coord + symop.displacement
end


canonize(arg::TranslationOperation) = iszero(arg.displacement) ? IdentityOperation() : arg

iscanonical(arg::TranslationOperation) = !iszero(arg.displacement)