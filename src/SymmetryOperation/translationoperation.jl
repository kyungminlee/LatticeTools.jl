export TranslationOperation

export inverse
export apply_symmetry
export canonize
export iscanonical

struct TranslationOperation{S<:Real} <:AbstractSymmetryOperation
    displacement::Vector{S}
    function TranslationOperation{S}(displacement::AbstractVector{S}) where {S<:Real}
        new{S}(displacement)
    end
    function TranslationOperation(displacement::AbstractVector{S}) where {S<:Real}
        new{S}(displacement)
    end
end

dimension(arg::TranslationOperation) = length(arg.displacement)

import Base.==
function (==)(lhs::TranslationOperation, rhs::TranslationOperation)
    return lhs.displacement == rhs.displacement
end

function (*)(lhs::TranslationOperation, rhs::TranslationOperation)
    return TranslationOperation(lhs.displacement .+ rhs.displacement)
end

function (*)(lhs::TranslationOperation, rhs::Real)
    return TranslationOperation(lhs.displacement .* rhs)
end

function (*)(lhs::Real, rhs::TranslationOperation)
    return TranslationOperation(lhs * rhs.displacement)
end

inverse(arg::TranslationOperation) = TranslationOperation(-arg.displacement)


function apply_symmetry(symop::TranslationOperation, coord::AbstractVector{<:Real})
    return coord + symop.displacement
end

function (symop::TranslationOperation)(coord::AbstractVector)
    return coord + symop.displacement
end


canonize(arg::TranslationOperation) = iszero(arg.displacement) ? IdentityOperation() : arg

iscanonical(arg::TranslationOperation) = !iszero(arg.displacement)